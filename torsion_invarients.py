import sympy as sp

Kpar, Kperp, beta = sp.symbols("Kpar Kperp beta", real=True)
H, a, W, Wr, r = sp.symbols("H a W Wr r", real=True)

gamma = 1 / sp.sqrt(1 - beta**2)

dim = 4
eta = sp.diag(-1, 1, 1, 1)

subs_general = {
    Kpar: H - Wr/a,
    Kperp: H - W/(a*r)
}

subs_static = {
    H: 0,
    a: 1
}


def zero_tensor():
    return sp.MutableDenseNDimArray.zeros(dim, dim, dim)


def simplify_value(value):
    if isinstance(value, list):
        return [sp.simplify(v) for v in value]
    return sp.simplify(value)


def substitute_value(value, subs):
    if isinstance(value, list):
        return [sp.simplify(v.subs(subs)) for v in value]
    return sp.simplify(value.subs(subs))


def build_initial_torsion():

    T = zero_tensor()

    eigenvalues = {
        1: Kpar,
        2: Kperp,
        3: Kperp
    }

    for i, Ki in eigenvalues.items():
        T[i, 0, i] = Ki
        T[i, i, 0] = -Ki

    return T


def radial_boost_matrix():

    return sp.Matrix([
        [gamma, -gamma*beta, 0, 0],
        [-gamma*beta, gamma, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])


def transform_torsion_covariantly(T):

    Lam = radial_boost_matrix()
    Lam_inv = Lam.inv()

    Tp = zero_tensor()

    for A in range(dim):
        for B in range(dim):
            for C in range(dim):

                expr = 0

                for D in range(dim):
                    for E in range(dim):
                        for F in range(dim):

                            expr += (
                                Lam[A, D]
                                * Lam_inv[E, B]
                                * Lam_inv[F, C]
                                * T[D, E, F]
                            )

                Tp[A, B, C] = sp.simplify(expr)

    return Tp


def lower_first_index(Tup):

    Tlow = zero_tensor()

    for A in range(dim):
        for B in range(dim):
            for C in range(dim):

                Tlow[A, B, C] = sp.simplify(
                    sum(
                        eta[A, D] * Tup[D, B, C]
                        for D in range(dim)
                    )
                )

    return Tlow


def raise_all_indices(Tlow):

    Tallup = zero_tensor()

    for A in range(dim):
        for B in range(dim):
            for C in range(dim):

                Tallup[A, B, C] = sp.simplify(
                    sum(
                        eta[A, D]
                        * eta[B, E]
                        * eta[C, F]
                        * Tlow[D, E, F]

                        for D in range(dim)
                        for E in range(dim)
                        for F in range(dim)
                    )
                )

    return Tallup


def torsion_vector(Tup):

    return [
        sp.simplify(
            sum(Tup[B, B, A] for B in range(dim))
        )
        for A in range(dim)
    ]


def vector_square(V):

    return sp.simplify(
        sum(
            eta[A, B] * V[A] * V[B]
            for A in range(dim)
            for B in range(dim)
        )
    )


def axial_vector(Tup):

    Tlow = lower_first_index(Tup)

    Avec = []

    for A in range(dim):

        expr = 0

        for B in range(dim):
            for C in range(dim):
                for D in range(dim):

                    expr += (
                        sp.LeviCivita(A, B, C, D)
                        * Tlow[B, C, D]
                    )

        Avec.append(sp.simplify(expr / 6))

    return Avec


def invariant_I1(Tup):

    Tlow = lower_first_index(Tup)
    Tallup = raise_all_indices(Tlow)

    return sp.simplify(

        sum(
            Tlow[A, B, C]
            * Tallup[A, B, C]

            for A in range(dim)
            for B in range(dim)
            for C in range(dim)
        )

    )


def invariant_I2(Tup):

    Tlow = lower_first_index(Tup)
    Tallup = raise_all_indices(Tlow)

    return sp.simplify(

        sum(
            Tlow[A, B, C]
            * Tallup[C, B, A]

            for A in range(dim)
            for B in range(dim)
            for C in range(dim)
        )

    )


def tegr_torsion_scalar(Tup):

    V = torsion_vector(Tup)

    return sp.simplify(

        sp.Rational(1, 4) * invariant_I1(Tup)
        + sp.Rational(1, 2) * invariant_I2(Tup)
        - vector_square(V)

    )


def all_invariants(Tup):

    V = torsion_vector(Tup)
    A = axial_vector(Tup)

    return {

        "I1 = T_abc T^abc":
            invariant_I1(Tup),

        "I2 = T_abc T^cba":
            invariant_I2(Tup),

        "T_a T^a":
            vector_square(V),

        "A_a A^a":
            vector_square(A),

        "TEGR torsion scalar T":
            tegr_torsion_scalar(Tup),

        "Torsion vector T_a":
            V,

        "Axial vector A^a":
            A,
    }


def print_block(title, data):

    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))

    for key, value in data.items():
        print(f"{key} = {simplify_value(value)}")


def print_substituted_block(title, data, subs):

    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))

    for key, value in data.items():
        print(f"{key} = {substitute_value(value, subs)}")


def print_differences(title, data1, data2):

    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))

    keys = [
        "I1 = T_abc T^abc",
        "I2 = T_abc T^cba",
        "T_a T^a",
        "A_a A^a",
        "TEGR torsion scalar T",
    ]

    for key in keys:

        print(
            f"{key}: "
            f"{sp.simplify(data2[key] - data1[key])}"
        )


T_initial = build_initial_torsion()

T_boosted = transform_torsion_covariantly(T_initial)

initial = all_invariants(T_initial)

boosted = all_invariants(T_boosted)


print_block(
    "Initial invariants: Kpar, Kperp form",
    initial
)

print_substituted_block(
    "Initial invariants: H, a, W, Wr, r form",
    initial,
    subs_general
)

print_substituted_block(
    "Initial invariants: static limit",
    initial,
    {**subs_general, **subs_static}
)

print_block(
    "Boosted invariants with spin connection",
    boosted
)

print_differences(
    "Boosted minus initial invariants",
    initial,
    boosted
)
