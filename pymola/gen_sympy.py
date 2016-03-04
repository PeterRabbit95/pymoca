#!/usr/bin/env python
"""
Modelica compiler.
"""
from __future__ import print_function
import sys
import antlr4
import antlr4.Parser
import argparse
import jinja2
from collections import OrderedDict
from .generated.ModelicaLexer import ModelicaLexer
from .generated.ModelicaParser import ModelicaParser
from .generated.ModelicaListener import ModelicaListener

#pylint: disable=invalid-name, no-self-use, missing-docstring, unused-variable, protected-access
#pylint: disable=too-many-public-methods

class SympyPrinter(ModelicaListener):

    #-------------------------------------------------------------------------
    # Setup
    #-------------------------------------------------------------------------

    def __init__(self, parser, trace):
        """
        Constructor
        """
        self._val_dict = OrderedDict()
        self.result = None
        self._parser = parser
        self._trace = trace
        self.indent = "            "

    @staticmethod
    def print_indented(ldr, s):
        if s is not None:
            for line in s.split('\n'):
                print(ldr, line)

    def setValue(self, ctx, val):
        """
        Sets tree values.
        """
        if ctx.depth() == 1:
            self.result = val
        else:
            self._val_dict[ctx] = val

    def getValue(self, ctx):
        """
        Gets tree values.
        """
        assert ctx is not None
        if ctx.depth() == 1:
            return self.result
        else:
            try:
                return self._val_dict[ctx]
            except KeyError:
                return None

    def enterEveryRule(self, ctx):
        if self._trace:
            ldr = " "*(ctx.depth()-1)
            rule = self._parser.ruleNames[ctx.getRuleIndex()]

            print(ldr, rule + "{")
            in_str = ctx.getText()
            if in_str > 0:
                print(ldr, "==============================")
                print(ldr, "INPUT\n")
                self.print_indented(ldr, ctx.getText())
                print(ldr, "==============================\n")

    def visitTerminal(self, node):
        pass

    def visitErrorNode(self, node):
        pass

    def exitEveryRule(self, ctx):
        rule = self._parser.ruleNames[ctx.getRuleIndex()]
        if self._trace:
            ldr = " "*ctx.depth()
            lines = None
            try:
                lines = self.getValue(ctx)
            except KeyError as e:
                pass

            if lines is not None:
                print(ldr, "==============================")
                print(ldr, "OUTPUT\n")
                self.print_indented(ldr, lines)
                print(ldr, "==============================\n")

            print(ldr + '} //' + rule + '\n')
        if self.getValue(ctx) is None:
            raise RuntimeError("no value set for {:s}:\ninput:\n{:s}".format(rule, ctx.getText()))

    #-------------------------------------------------------------------------
    # Top Level
    #-------------------------------------------------------------------------

    def exitStored_definition(self, ctx):
        composition = self.getValue(ctx.class_definition()[0]\
                .class_specifier().composition())
        self.result = """
#############################################################################
# Automatically generated by pymola

from __future__ import print_function, division
import sympy
assert sympy.__version__ >= '0.7.6.1'
import sympy.physics.mechanics as mech
sympy.init_printing()
try:
    mech.init_vprinting()
except AttributeError:
    mech.mechanics_printing()
import scipy.integrate
import pylab as pl
from collections import OrderedDict

#pylint: disable=too-few-public-methods, too-many-locals, invalid-name, no-member

class Model(object):
    \"\"\"
    Modelica Model.
    \"\"\"

    def __init__(self):
        \"\"\"
        Constructor.
        \"\"\"

        self.t = sympy.symbols('t')

        {composition:s}

        self.x = sympy.Matrix(self.x)
        self.x_dot = self.x.diff(self.t)

        self.sol = sympy.solve(self.eqs, self.x_dot)

        self.f = sympy.Matrix([self.sol[xi] for xi in self.x_dot])
        print('x:', self.x)
        print('f:', self.f)

        self.p_vect = [locals()[key] for key in self.p_dict.keys()]

        print('p:', self.p_vect)

        self.f_lam = sympy.lambdify((self.t, self.x, self.p_vect), self.f)

    def get_p0(self):
        return [self.p_dict[key] for key in
            sorted(self.p_dict.keys())]

    def get_x0(self):
        return [self.x0_dict[key] for key in
            sorted(self.x0_dict.keys())]

    def simulate(self, tf=30, dt=0.001, show=False):
        \"\"\"
        Simulation function.
        \"\"\"

        p0 = self.get_p0()
        x0 = self.get_x0()

        print('p0', p0)
        print('x0', x0)

        sim = scipy.integrate.ode(self.f_lam)
        sim.set_initial_value(x0, 0)
        sim.set_f_params(p0)

        data = {{
            'x': [],
            't': [],
        }}

        while  sim.t < tf:
            sim.integrate(sim.t + dt)
            t = sim.t

            # TODO replace hardcoded when statement
            # below
            #velocity = sim.y[0]
            #height = sim.y[1]
            #c = self.p0[0]
            #if velocity < 0 and height < 0:
            #    velocity = -c*velocity
            #    height = 0
            #    sim.set_initial_value([velocity, height], t)

            # data['x'] += [[velocity, height]]
            data['x'] += [sim.y]
            data['t'] += [t]

        pl.plot(data['t'], data['x'])
        if show:
            pl.show()

        return data

if __name__ == "__main__":
    model = Model()
    model.simulate()

#############################################################################
""".format(**locals())

    def exitComposition(self, ctx):
        decl = self.getValue(ctx.element_list()[0])
        equations = self.getValue(ctx.equation_section()[0])
        data = locals()
        data['walker'] = self
        data.pop('self')
        self.setValue(ctx, jinja2.Template("""
        {%- set params=decl.parameters -%}
        {%- set dsyms=decl.dynamic_symbols -%}
        {%- set start_values=decl.start_values -%}

        {% if params|length > 0 %}
        # symbols
        {{params|join(', ') }} = \\
            sympy.symbols('{{params|join(', ')}}')
        {% endif %}

        # dynamic symbols
        {{dsyms|join(', ')}} = \\
            mech.dynamicsymbols('{{dsyms|join(', ')}}')

        # parameters
        self.p_dict = OrderedDict({
        {%- for key in params.keys() %}
            '{{key}}': {{params[key]}},
        {%- endfor %}
        })

        # initial sate
        self.x0_dict = OrderedDict({
        {%- for key in dsyms.keys() %}
            '{{key}}': {{dsyms[key].start}},
        {%- endfor %}
        })

        # state space
        self.x = sympy.Matrix([
            {{dsyms|join(', ')}}
        ])
{{equations}}
""").render(**data))

    def exitElement_list(self, ctx):
        d = self.getValue(ctx.element()[0])
        for i in range(1, len(ctx.element())):
            element = ctx.element()[i]
            for key in d.keys():
                d[key].update(self.getValue(element)[key])
        self.setValue(ctx, d)
        print("dict", d)

    def exitElement(self, ctx):
        self.setValue(ctx, self.getValue(ctx.component_clause()))

    def exitEquation_section(self, ctx):
        str_eq = "self.eqs = ["
        str_when = ""
        for eq in ctx.equation():
            data = self.getValue(eq)
            if len(data) <= 0:
                continue
            print("eq type:", type(eq.equation_options()))
            if isinstance(eq.equation_options(), ModelicaParser.Equation_simpleContext):
                str_eq += "\n{self.indent:s}{data:s},".format(**locals())
                print("SIMPLE")
            elif isinstance(eq.equation_options(), ModelicaParser.Equation_whenContext):
                str_when += "\n{self.indent:s}{data:s}".format(**locals())
            else:
                raise IOError('equation type not supported yet')
        str_eq += "\n{self.indent:s}]\n".format(**locals())
        str_when += ""
        self.setValue(ctx, """
        # equations
        {str_eq:s}

        # when equations
{str_when:s}
""".format(**locals()))

    def exitComponent_clause(self, ctx):
        s = ""
        parameters = OrderedDict()
        dynamic_symbols = OrderedDict()

        if ctx.type_prefix().getText() == 'parameter':
            # store all parameters
            for comp in ctx.component_list().component_declaration():
                name = comp.declaration().IDENT().getText()
                val = comp.declaration().modification().expression().getText()
                parameters[name] = val
        else:
            # store all variables
            for comp in ctx.component_list().component_declaration():
                name = comp.declaration().IDENT().getText()
                dynamic_symbols[name] = {'start': 0}
                try:
                    mod = comp.declaration().modification().class_modification()
                    arg = mod.argument_list().argument()[0]
                    emod = arg.element_modification_or_replaceable().element_modification()
                    if emod.name().getText() == 'start':
                        val = emod.modification().expression().getText()
                        dynamic_symbols[name]['start'] = val
                except AttributeError:
                    pass

        self.setValue(ctx, {
            'parameters': parameters,
            'dynamic_symbols': dynamic_symbols,
            })

    def exitClass_prefixes(self, ctx):
        # don't care about prefixes for now
        self.setValue(ctx, '')

    def exitName(self, ctx):
        self.setValue(ctx, ctx.getText())

    def exitType_specifier(self, ctx):
        self.setValue(ctx, self.getValue(ctx.name()))

    def exitModification_class(self, ctx):
        mod = self.getValue(ctx.class_modification())
        if ctx.expression() is not None:
            expr = self.getValue(ctx.expression())
            self.setValue(ctx, "{:s} = {:s}".format(mod, expr))
        else:
            self.setValue(ctx, "{:s}".format(mod))

    def exitModification_assignment(self, ctx):
        expr = self.getValue(ctx.expression())
        self.setValue(ctx, "= {:s}".format(expr))

    def exitModification_assignment2(self, ctx):
        # TODO how to handle := operator?, is it just assignment
        expr = self.getValue(ctx.expression())
        self.setValue(ctx, "= {:s}".format(expr))

    def exitDeclaration(self, ctx):
        # set above based on component type
        self.setValue(ctx, "")

    def exitComment(self, ctx):
        self.setValue(ctx, "# {:s}".format(ctx.getText()))

    def exitString_comment(self, ctx):
        self.setValue(ctx, "# {:s}".format(ctx.getText()))

    def exitComponent_declaration(self, ctx):
        # set above based on component type
        self.setValue(ctx, "")

    def exitComponent_list(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitElement_modification(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitElement_modification_or_replaceable(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitType_prefix(self, ctx):
        self.setValue(ctx, "")

    def exitArgument(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitArgument_list(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitClass_modification(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitComponent_reference(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitFunction_arguments(self, ctx):
        if ctx.function_argument() is not None:
            s = ""
            args = ctx.function_argument()
            n = len(args)
            for i in range(n):
                s += "{:s}".format(self.getValue(args[i]))
                if i < n-1:
                    s += ", "
            self.setValue(ctx, s)
        elif ctx.named_arguments() is not None:
            args = self.getValue(ctx.named_arguments())
            self.setValue(ctx, args)

    def exitArguement_function(self, ctx):
        name = self.getValue(ctx.name())
        if ctx.named_arguments() is not None:
            named_arguments = self.getValue(ctx.named_arguements())
        else:
            named_arguments = ""
        self.setValue(ctx, "{name:s}({named_arguments:s})")

    def exitArgument_expression(self, ctx):
        expr = self.getValue(ctx.expression())
        self.setValue(ctx, expr)

    def exitFunction_call_args(self, ctx):
        args = ctx.function_arguments()
        if  args is not None:
            self.setValue(ctx, self.getValue(args))
        else:
            self.setValue(ctx, "")

    def exitClass_specifier(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitClass_definition(self, ctx):
        # TODO
        self.setValue(ctx, "")

    def exitClass_stored_definition(self, ctx):
        # TODO
        self.setValue(ctx, "")

    #-------------------------------------------------------------------------
    # Equation
    #-------------------------------------------------------------------------

    def exitEquation(self, ctx):
        self.setValue(ctx, self.getValue(ctx.equation_options()))

    def exitEquation_simple(self, ctx):
        self.setValue(
            ctx,
            "{:s} - {:s}".format(
                self.getValue(ctx.simple_expression()),
                self.getValue(ctx.expression())))

    def exitEquation_if(self, ctx):
        raise NotImplementedError("")

    def exitEquation_for(self, ctx):
        raise NotImplementedError("")

    def exitEquation_connect_clause(self, ctx):
        raise NotImplementedError("")

    def exitEquation_function(self, ctx):
        self.setValue(ctx, "{name:s}({args:})".format(**{
            'name': self.getValue(ctx.name()),
            'args': self.getValue(ctx.function_call_args())
            }))

    def exitWhen_equation(self, ctx):
        exprs = ctx.expression()
        eqns = ctx.equation()
        t = jinja2.Template("""
if {{ walker.getValue(exprs[0]) -}}:
    {{ walker.getValue(eqns[0]) }}
{% for i in range(1, exprs|length) %}
elif {{ walker.getValue(exprs[i]) }}:
    {{ walker.getValue(eqns[i]) }}
{% endfor %}
""")
        # self.setValue(ctx, t.render({
            # 'walker': self,
            # 'exprs': exprs,
            # 'eqns': eqns,
            # }))
        self.setValue(ctx, '')

    def exitEquation_when(self, ctx):
        self.setValue(ctx, self.getValue(ctx.when_equation()))

    #-------------------------------------------------------------------------
    # Expression
    #-------------------------------------------------------------------------

    def exitExpression_simple(self, ctx):
        self.setValue(ctx, self.getValue(ctx.simple_expression()))

    def exitSimple_expression(self, ctx):
        # TODO, can have more than one expr
        self.setValue(ctx, self.getValue(ctx.expr()[0]))

    def exitExpr_primary(self, ctx):
        self.setValue(ctx, self.getValue(ctx.primary()))

    def exitExpr_neg(self, ctx):
        self.setValue(ctx, '-({:s})'.format(self.getValue(ctx.expr())))

    def exitExpr_exp(self, ctx):
        self.setValue(ctx, '(({:s})**({:s}))'.format(
            self.getValue(ctx.primary()[0]),
            self.getValue(ctx.primary()[1])))

    def exitExpr_rel(self, ctx):
        self.setValue(ctx, '({:s} {:s} {:s})'.format(
            self.getValue(ctx.expr()[0]),
            ctx.op.text,
            self.getValue(ctx.expr()[1])))

    def exitExpr_mul(self, ctx):
        self.setValue(ctx, '({:s} {:s} {:s})'.format(
            self.getValue(ctx.expr()[0]),
            ctx.op.text,
            self.getValue(ctx.expr()[1])))

    def exitExpr_add(self, ctx):
        self.setValue(ctx, '({:s} {:s} {:s})'.format(
            self.getValue(ctx.expr()[0]),
            ctx.op.text,
            self.getValue(ctx.expr()[1])))

    #-------------------------------------------------------------------------
    # Primary
    #-------------------------------------------------------------------------

    def exitPrimary_unsigned_number(self, ctx):
        self.setValue(ctx, ctx.getText())

    def exitPrimary_component_reference(self, ctx):
        self.setValue(ctx, ctx.getText())

    def exitPrimary_derivative(self, ctx):
        name = ctx.function_call_args().function_arguments().function_argument()[0].getText()
        self.setValue(ctx, '{:s}.diff(self.t)'.format(name))

    def exitPrimary_string(self, ctx):
        self.setValue(ctx, ctx.getText())

    def exitPrimary_false(self, ctx):
        self.setValue(ctx, ctx.getText())


def generate(modelica_model, trace=False):
    "The modelica model"
    input_stream = antlr4.InputStream(modelica_model)
    lexer = ModelicaLexer(input_stream)
    stream = antlr4.CommonTokenStream(lexer)
    parser = ModelicaParser(stream)
    tree = parser.stored_definition()
    # print(tree.toStringTree(recog=parser))
    sympyPrinter = SympyPrinter(parser, trace)
    walker = antlr4.ParseTreeWalker()
    walker.walk(sympyPrinter, tree)
    return sympyPrinter.result


def main(argv):
    #pylint: disable=unused-argument
    "The main function"
    parser = argparse.ArgumentParser()
    parser.add_argument('src')
    parser.add_argument('out')
    parser.add_argument('-t', '--trace', action='store_true')
    parser.set_defaults(trace=False)
    args = parser.parse_args(argv)
    with open(args.src, 'r') as f:
        modelica_model = f.read()
    sympy_model = generate(modelica_model, trace=args.trace)

    with open(args.out, 'w') as f:
        f.write(sympy_model)

if __name__ == '__main__':
    main(sys.argv[1:])

# vim: set et ft=python fenc=utf-8 ff=unix sts=0 sw=4 ts=4 :
