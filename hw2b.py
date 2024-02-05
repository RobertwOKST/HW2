from math import cos, pi


def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    This function implements the Secant method to find the root of an equation.
    :param fcn: the function for which we want to find the root
    :param x0: x value in the neighborhood of the root (or guess 1)
    :param x1: another x value in the neighborhood of the root (or guess x0+1)
    :param maxiter: exit if the number of iterations (new x values) equals this number
    :param xtol: exit if the |xnewest - xprevious| < xtol
    :return: tuple with: (the final estimate of the root (most recent value of x), number of iterations)
    """
    for i in range(maxiter):
        fx0 = fcn(x0)
        fx1 = fcn(x1)

        if abs(fx1 - fx0) < xtol:
            return x1, i

        x_new = x1 - (fx1 * (x1 - x0)) / (fx1 - fx0)

        x0, x1 = x1, x_new

    return x1, maxiter


def fn1(x):
    return x - 3 * cos(x)


def fn2(x):
    return cos(2 * x) * x ** 3


def main():
    """
    fn1: x - 3cos(x) = 0; with x0=1, x1=2, maxiter=5, and xtol=1e-4
    fn2: cos(2x)*x**3; with x0=1, x1=2, maxiter=15, and xtol=1e-8
    fn3: cos(2x)*x**3; with x0=1, x1=2, maxiter=3, and xtol=1e-8

    I observe that for functions 2 and 3, the answer should be pi/2 or about 1.57
    :return: nothing, just print results
    """
    r1 = Secant(fn1, 1, 2, 5, 1e-4)
    r2 = Secant(fn2, 1, 2, 15, 1e-8)
    r3 = Secant(fn2, 1, 2, 3, 1e-8)

    print("Root of fn1 = {root:.4f}, after {iter:d} iterations".format(root=r1[0], iter=r1[1]))
    print("Root of fn2 = {root:.4f}, after {iter:d} iterations".format(root=r2[0], iter=r2[1]))
    print("Root of fn3 = {root:.4f}, after {iter:d} iterations".format(root=r3[0], iter=r3[1]))


if __name__ == "__main__":
    main()