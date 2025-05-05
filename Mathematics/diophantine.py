from sympy import symbols, Eq, diophantine
from itertools import product
from _helpers import pprint_format_as_set


def solve_linear_diophantine_nonnegative(
    coeffs: list[int],
    constant: int,
    var_names: list[str] = None,
    verbose: bool = True
) -> list[tuple[int, ...]]:
    """
    Solves a linear Diophantine equation of the form:
        a1*x1 + a2*x2 + ... + an*xn = c
    With added restriction for results to be nonnegative.

    Args:
        coeffs (list[int]):
            List of coefficients [a1, a2, ..., an].
        constant (int):
            Target value to solve.
        var_names (list[str]):
            Optional list of variable names ['x', 'y', 'z', ...]. Defaults to None.

    Returns:
        list[tuple[int, ...]]: List of tuples containing nonnegative integer solutions.
    """

    # Validate coefficients
    for coeff in coeffs:
        assert coeff != 0, 'Each coefficient must be nonzero.'

    # Get number of coefficients
    n = len(coeffs)

    # Resolve variable names
    if var_names is None:
        var_names = [f'x{i+1}' for i in range(n)]

    # Deduplicate identical coefficients, prioritizing first occurrence
    seen = set()
    simplified_coeffs, simplified_vars = [], []
    for coeff, var_name in zip(coeffs, var_names):
        if coeff not in seen:
            seen.add(coeff)
            simplified_coeffs.append(coeff)
            simplified_vars.append(var_name)
    if len(simplified_vars) < len(var_names) and verbose:
        print('[INFO] Simplified equation due to duplicate coefficients')

    # Instantiante Diophantine equation
    var_symbols = symbols(simplified_vars)
    lhs = sum(coeff * var for coeff, var in zip(simplified_coeffs, var_symbols))
    eq = Eq(lhs, constant)

    # Obtain analytical parametric solution
    sol_set = diophantine(eq)
    if verbose:
        print(f'[INFO] Analytical solution is {sol_set}')

    # Initialize solution set then iterate
    valid_solutions = set()
    for sol in sol_set:
        # Get all parameters used across all expressions
        all_params = sorted(set().union(*[expr.free_symbols for expr in sol]), key=lambda s: s.name)
        
        if not all_params:
            # Handle finite valid solution
            if all(expr.is_integer and expr >= 0 for expr in sol):
                valid_solutions.add(tuple(int(expr) for expr in sol))
        else:
            # Get minimum magnitude for parameter values
            min_coeff = min(abs(coeff) for coeff in simplified_coeffs)
            # Calculate effective bound (sum of minimum term contributions <= constant)
            param_bound = 3 * constant // min_coeff  # TODO: Improve parameter bound calculation logic
            # Iterate to append valid solutions
            for values in product(range(-param_bound, param_bound + 1), repeat=len(all_params)):
                subs = dict(zip(all_params, values))
                eval_sol = [expr.subs(subs) for expr in sol]
                if all(expr.is_integer and expr >= 0 for expr in eval_sol):
                    valid_solutions.add(tuple(int(expr) for expr in eval_sol))
    
    # Sort solutions
    result = sorted(valid_solutions)

    # Return variables and result
    return result, simplified_vars


if __name__ == '__main__':
    # Example Usage:
    coefficients = [11, 8, 9]
    variables = ['x', 'y', 'z']
    target = 96

    # Solve 11x + 8y + 9z = 96
    solutions, solution_vars = solve_linear_diophantine_nonnegative(coefficients, target, variables)
    print('[RESULT]')
    pprint_format_as_set(solutions, solution_vars)
