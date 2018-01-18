import casadi as ca
import numpy as np

import os

import unittest
import pymola.backends.casadi.generator as gen_casadi
from pymola.backends.casadi.alias_relation import AliasRelation
from pymola.backends.casadi.model import Model, Variable
from pymola.backends.casadi.api import transfer_model, CachedModel
from pymola import parser, ast

TEST_DIR = os.path.dirname(os.path.realpath(__file__))

class GenCasadiTest(unittest.TestCase):
    def assert_model_equivalent(self, A, B):
        def sstr(a): return set([str(e) for e in a])

        for l in ["states", "der_states", "inputs", "outputs", "constants", "parameters"]:
            self.assertEqual(sstr(getattr(A, l)), sstr(getattr(B, l)))

    def assert_model_equivalent_numeric(self, A, B, tol=1e-9):
        self.assertEqual(len(A.states), len(B.states))
        self.assertEqual(len(A.der_states), len(B.der_states))
        self.assertEqual(len(A.inputs), len(B.inputs))
        self.assertEqual(len(A.outputs), len(B.outputs))
        self.assertEqual(len(A.constants), len(B.constants))
        self.assertEqual(len(A.parameters), len(B.parameters))

        if not isinstance(A, CachedModel) and not isinstance(B, CachedModel):
            self.assertEqual(len(A.equations), len(B.equations))
            self.assertEqual(len(A.initial_equations), len(B.initial_equations))

        for f_name in ['dae_residual', 'initial_residual', 'variable_metadata']:
            this = getattr(A, f_name + '_function')
            that = getattr(B, f_name + '_function')

            np.random.seed(0)

            args_in = []
            for i in range(this.n_in()):
                sp = this.sparsity_in(0)
                r = ca.DM(sp, np.random.random(sp.nnz()))
                args_in.append(r)

            this_out = this.call(args_in)
            that_out = that.call(args_in)

            # N.B. Here we require that the order of the equations in the two models is identical.
            for i, (a, b) in enumerate(zip(this_out, that_out)):
                for j in range(a.size1()):
                    for k in range(a.size2()):
                        if a[j, k].is_regular() or b[j, k].is_regular():
                            test = float(ca.norm_2(ca.vec(a[j, k] - b[j, k]))) <= tol
                            if not test:
                                print(j)
                                print(k)
                                print(a[j,k])
                                print(b[j,k])
                                print(f_name)
                            self.assertTrue(test)

        return True


    def test_simple_array(self):
        with open(os.path.join(TEST_DIR, 'SimpleArray.mo'), 'r') as f:
            txt = f.read()
        ast_tree = parser.parse(txt)
        casadi_model = gen_casadi.generate(ast_tree, 'SimpleArray')

        casadi_model.simplify({})

        print(casadi_model)
        ref_model = Model()

        a_1 = ca.MX.sym("a_1")
        a_2 = ca.MX.sym("a_2")
        a_3 = ca.MX.sym("a_3")
        b_1 = ca.MX.sym("b_1")
        b_2 = ca.MX.sym("b_2")
        b_3 = ca.MX.sym("b_3")
        b_4 = ca.MX.sym("b_4")
        c_1 = ca.MX.sym("c_1")
        c_2 = ca.MX.sym("c_2")
        c_3 = ca.MX.sym("c_3")

        ref_model.alg_states = list(map(Variable, [a_1, a_2, a_3, c_1, c_2, c_3]))

        for i in range(3, 6):
            ref_model.alg_states[i].min = 0.0

        # No parameters
        ref_model.constants = list(map(Variable, [b_1, b_2, b_3, b_4]))
        constant_values = [2.7, 3.7, 4.7, 5.7]
        for const, val in zip(ref_model.constants, constant_values):
            const.value = val

        ref_model.equations = [c_1 - (a_1 + b_1),
                               c_2 - (a_2 + b_2),
                               c_3 - (a_3 + b_3),
                               a_1 - 1,
                               a_2 - 2,
                               a_3 - 3]

        if not self.assert_model_equivalent_numeric(ref_model, casadi_model):
            raise Exception("Failed test")


c = GenCasadiTest()
c.test_simple_array()