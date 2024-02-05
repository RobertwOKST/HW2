from copy import deepcopy

def GaussSeidel(Aaug, x, Niter=15):
    """
    This implements the Gauss-Seidel method for solving a system of equations.
    :param Aaug: The augmented matrix from Ax=b -> [A|b]
    :param x: An initial guess for the x vector. if A is nxn, x is nx1
    :param Niter: Number of iterations to run the GS method
    :return: the solution vector x
    """

    A = [row[:-1] for row in Aaug]  # Extract coefficient matrix A
    b = [row[-1] for row in Aaug]   # Extract the right-hand side vector b

    n = len(b)

    # Make the matrix diagonal dominant
    A = MakeDiagDom(A)

    for _ in range(Niter):
        for i in range(n):
            sum_ax = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - sum_ax) / A[i][i]

    return x

def MakeDiagDom(A):
    """
    This function reorders the rows of matrix A to put the largest absolute values along the diagonal.
    :param A: The matrix to sort
    :return: The sorted matrix
    """

    n = len(A)
    B = deepcopy(A)

    for i in range(n):
        max_index = max(range(i, n), key=lambda j: abs(B[j][i]))
        if i != max_index:
            B[i], B[max_index] = B[max_index], B[i]

    return B

def main():
    """This portion sets the matrices for equation 1 and equation 2
    to solve for the x values"""
    #Equation 1
    Aaug1 = [[3, 1, -1, 2],
             [1, 4, 1, 12],
             [2, 1, 2, 10]]

    x1 = [0.0] * 3  # Initial guess for example 1
    solution1 = GaussSeidel(Aaug1, x1)
    print("Solution for Equation 1:")
    print("x1 = {A:0.3f}".format(A=solution1[0]))
    print("x2 = {A:0.3f}".format(A=solution1[1]))
    print("x3 = {A:0.3f}".format(A=solution1[2]))

    # Equation 2
    Aaug2 = [[1, -10, 2, 4, 2],
             [3, 1, 4, 12, 12],
             [9, 2, 3, 4, 21],
             [-1, 2, 7, 3, 37]]

    x2 = [0.0] * 4  # Initial guess for example 2
    solution2 = GaussSeidel(Aaug2, x2)
    print("Solution for Equation 2:")
    print("x1 = {A:0.3f}".format(A=solution2[0]))
    print("x2 = {A:0.3f}".format(A=solution2[1]))
    print("x3 = {A:0.3f}".format(A=solution2[2]))
    print("x4 = {A:0.3f}".format(A=solution2[3]))

if __name__ == "__main__":
    main()