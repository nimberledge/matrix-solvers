"""Helper file to generate matrices for benchmark tests"""
import numpy as np
import scipy.linalg as la
import numpy.linalg as nla
from sklearn.datasets import make_spd_matrix
from sklearn.datasets import make_sparse_spd_matrix


# Generates an spd matrix of size NxN
# For dense cholesky
def generate_spd_matrix(N):
    X = make_spd_matrix(N, random_state=0)
    return X


# Write a matrix in the format we laid out
def write_to_file(M, filename):
    with open(filename, "w+") as in_file:
        first_line = str(M.shape[0]) + "\t" + str(M.shape[1]) + "\n"
        second_line = "\t".join([str(i) for i in M.flatten()]) + "\n"
        in_file.writelines([first_line, second_line])
    print("Data written to: " + filename)


# Use this for dense Jacobi
def generate_jac_matrix(N):
    # Generate N random eigenvalues < 1
    eigvals = np.random.random(N)

    S = np.diag(eigvals)
    q, _ = la.qr(np.random.rand(N, N))
    semidef = q.T @ S @ q
    for i in range(N):
        # Make diagonally dominant
        semidef[i, i] = sum(abs(semidef[i, :])) + 0.0001

    return semidef


# Use this for the sparse jacobi
def generate_sparse_spd(N):
    X = make_sparse_spd_matrix(dim=N, random_state=0)
    for i in range(N):
        # Make diagonally dominant
        X[i, i] = sum(abs(X[i, :])) + 0.0001
    return X


# Write in our laid out sparse format
def write_sparse_to_file(M, filename):
    nnz = 0
    zero_tol = 1e-10
    pts = []
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if abs(M[i][j]) > zero_tol:
                nnz += 1
                pts.append({'row': i, 'col': j, 'val': M[i][j]})
    first_line = str(M.shape[0]) + "\t" + str(M.shape[1]) + "\t" + str(nnz) + "\n"  # noqa
    lines = [first_line]
    for pt in pts:
        line = str(pt['row']) + "\t" + str(pt['col']) + "\t" + str(pt["val"]) + "\n"  # noqa
        lines.append(line)
    with open(filename, "w+") as out_file:
        out_file.writelines(lines)
    print("Data written to: " + filename)


def main():
    # print("Hello there, really sorry this file exists")
    # np.random.seed(5)
    sizes = [40]
    for N in sizes:
        M = generate_sparse_spd(N)
        # Condition number was a half-decent indicator
        # of if I'd messed up
        print(nla.cond(M))
        filename = "sparsetest{}by{}.txt".format(N, N)
        write_sparse_to_file(M, filename)
        # write_to_file(M, filename)


if __name__ == '__main__':
    main()
