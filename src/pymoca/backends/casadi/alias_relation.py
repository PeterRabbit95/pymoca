from collections import OrderedDict
import casadi as ca


OPERATION = {
    ca.OP_MUL: lambda x1, x2: x1 * x2,
    ca.OP_DIV: lambda x1, x2: x1 / x2,
    ca.OP_NEG: lambda x: -x,
    ca.OP_PARAMETER: lambda x: x,
    ca.OP_TWICE: lambda x: 2*x,
    'OP_TWICE_INVERSE': lambda x: x/2,
    ca.OP_SIN: lambda x: ca.sin(x),
    ca.OP_COS: lambda x: ca.cos(x),
    ca.OP_ASIN: lambda x: ca.asin(x),
    ca.OP_ACOS: lambda x: ca.acos(x)
}

INVERSE_OPERATOR = {
    ca.OP_MUL: ca.OP_DIV,
    ca.OP_DIV: ca.OP_MUL,
    ca.OP_NEG: ca.OP_PARAMETER,
    ca.OP_PARAMETER: ca.OP_NEG,
    ca.OP_TWICE: 'OP_TWICE_INVERSE',
    ca.OP_COS: ca.OP_ACOS,
    ca.OP_SIN: ca.OP_ASIN,
}


def casadi_operation(operands, operator, inverse = False):
    if inverse:
        operator = INVERSE_OPERATOR[operator]

    return OPERATION[operator](*operands)


class Alias:
    def __init__(self, alias, canonical, direct, inverse):
        self.alias = alias
        self.canonical = canonical
        self.direct = direct
        self.inverse = inverse
        self.direct_call = ca.Function('Direct', [self.canonical], [direct])
        self.inverse_call = ca.Function('Inverse', [self.alias], [inverse])

    def __repr__(self):
        return 'Can: {}\nAlias: {}\nDirect: {}\nInverse: {}'.format(self.canonical, self.alias, self.direct, self.inverse)


# Code snippet from RTC-Tools, copyright Stichting Deltares, originally under the terms of the GPL
# version 3.  Relicensed with permission.
class AliasRelation:
    def __init__(self):
        self._aliases = []

    def remove_canonical_variable(self, var):
        self._aliases[:] = [alias for alias in self._aliases if not ca.depends_on(alias.canonical, var)]

    @property
    def canonical_variables(self):
        canonical_variables = set()
        for alias in self._aliases:
            canonical_variables.add(alias.canonical)
        return canonical_variables

    def __iter__(self):
        for canonical_variable in self.canonical_variables:
            aliases = [alias for alias in self._aliases if ca.depends_on(alias.canonical, canonical_variable)]
            yield canonical_variable, aliases

    def alias_from_equation(self, eq: ca.MX, canonical, alias):
        try:
            # Separate equation elements
            if ca.depends_on(eq.dep(0), alias) and ca.depends_on(eq.dep(1), canonical):
                left0 = eq.dep(0)
                left1 = eq.dep(1)
            elif ca.depends_on(eq.dep(1), alias) and ca.depends_on(eq.dep(0), canonical):
                left0 = eq.dep(1)
                left1 = eq.dep(0)
            else:
                return False

            # Invert equations
            direct = self._explicit_equation_2_dep(left0, left1, eq.op())
            inverse = self._explicit_equation_2_dep(left1, left0, eq.op())

            # Build alias
            self._aliases.append(Alias(alias, canonical, direct, inverse))
            return True
        except:
            return False

    @staticmethod
    def _explicit_equation_2_dep(left0, left1, main_op):
        # Handle main operation (move left term 1 to the right changing sign)
        if main_op == ca.OP_SUB:
            right = left1
        elif main_op == ca.OP_ADD:
            right = -left1
        else:
            raise NotImplementedError('Cannot handle main operations which are not addition or subtractions')

        # Change sign if the term 0 is negative
        if left0.n_dep() == 1 and left0.op() == ca.OP_NEG:
            left0 = -left0
            right = -right

        if left0.n_dep() == 0:  # If left is single return right - all operations have been performed
            # Only considered case is parameter now
            assert left0.op() == ca.OP_PARAMETER
            return right

        elif left0.n_dep() == 1:  # If there is one more operation apply it
            if left0.dep(0).op() == ca.OP_PARAMETER or left0.dep(0).dep(0) == ca.OP_PARAMETER:
                return OPERATION[INVERSE_OPERATOR[left0.op()]](right)
            else:
                raise NotImplementedError

        elif left0.n_dep() == 2:

            # Isolate numeric part
            if left0.dep(0).is_constant() and left0.dep(1).is_symbolic():
                if left0.op() == ca.OP_DIV:
                    right = 1 / right
                num = left0.dep(0)

            elif left0.dep(1).is_constant() and left0.dep(0).is_symbolic():
                num = left0.dep(1)

            else:
                raise NotImplementedError

            return casadi_operation([right, num], left0.op(), True)

# # test
# x = ca.MX.sym('x')
# y = ca.MX.sym('y')
#
# alias = AliasRelation()
# alias.alias_from_equation(5*x - 9*y, x, y)
# alias.alias_from_equation(2*x - 9*y, x, y)
# alias.alias_from_equation(x - 9*y, x, y)
# alias.alias_from_equation(x - y, x, y)
# alias.alias_from_equation(-x - y, x, y)
#
# for al, rel in alias:
#     print('{}: {}'.format(al, rel))
