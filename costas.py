import numpy as np
from numba import njit
import matplotlib.pyplot as plt


# Optimized utility functions
@njit
def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for x in range(3, int(n**0.5) + 1, 2):
        if n % x == 0:
            return False
    return True


@njit
def powmod(base, exp, mod):
    result = 1
    while exp > 0:
        if exp & 1:
            result = (result * base) % mod
        base = (base * base) % mod
        exp >>= 1
    return result


@njit
def find_primitive_root(p):
    phi = p - 1
    factors = []
    temp = phi
    for f in range(2, int(temp**0.5) + 1):
        if temp % f == 0:
            factors.append(f)
            while temp % f == 0:
                temp //= f
    if temp > 1:
        factors.append(temp)

    for g in range(2, p):
        for factor in factors:
            if powmod(g, phi // factor, p) == 1:
                break
        else:
            return g
    return -1


@njit
def welch_costas(n):
    p = n + 1
    if not is_prime(p):
        return None

    g = find_primitive_root(p)
    if g < 0:
        return None

    perm = np.empty(n, dtype=np.int64)
    cur = 1
    for i in range(n):
        cur = (cur * g) % p
        perm[i] = cur - 1
    return perm


# Optimized backtracking search
@njit
def _search_costas(permutation, used_col, used_diffs, row, N):
    if row == N:
        return True

    for col in range(N):
        if not used_col[col]:
            valid = True
            for rprev in range(row):
                rdiff = row - rprev
                cdiff = col - permutation[rprev] + N - 1
                if used_diffs[rdiff, cdiff]:
                    valid = False
                    break

            if valid:
                permutation[row] = col
                used_col[col] = True
                for rprev in range(row):
                    rdiff = row - rprev
                    cdiff = col - permutation[rprev] + N - 1
                    used_diffs[rdiff, cdiff] = True

                if _search_costas(permutation, used_col, used_diffs, row + 1, N):
                    return True

                used_col[col] = False
                for rprev in range(row):
                    rdiff = row - rprev
                    cdiff = col - permutation[rprev] + N - 1
                    used_diffs[rdiff, cdiff] = False

    return False


@njit
def costas_array_backtracking(N):
    if N == 1:
        return np.array([0], dtype=np.int64)

    permutation = np.full(N, -1, dtype=np.int64)
    used_col = np.zeros(N, dtype=np.bool)
    used_diffs = np.zeros((N, 2 * N - 1), dtype=np.bool)

    if _search_costas(permutation, used_col, used_diffs, 0, N):
        return permutation
    return None


@njit
def generate_costas_array_njit(N):
    if N < 1:
        return None

    permutation = welch_costas(N)
    if permutation is not None:
        return permutation

    return costas_array_backtracking(N)


# Combined generator
def generate_costas_array(N):
    result = generate_costas_array_njit(N)
    if result is not None:
        (
            print("Generated using Welch method.")
            if welch_costas(N) is not None
            else print("Falling back to backtracking search...")
        )
    return result


# Visualization function
def visualize_costas_array(permutation):
    if permutation is None:
        print("No Costas array to display.")
        return

    N = len(permutation)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(np.arange(N) + 0.5, permutation + 0.5, color="blue", s=100 / N)
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"Costas Array (N={N})")
    plt.gca().invert_yaxis()
    plt.show()
