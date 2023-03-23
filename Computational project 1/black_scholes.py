"""
This file contains the numerical schemes for solving the Black-Scholes equation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, animation
from mpl_toolkits.mplot3d import axes3d


def apply_black_scholes(strike_price, put: bool, boundary, time, space, interest_rate, n_t, n_s, sigma, black_scholes):
    return black_scholes(strike_price, put, boundary, time, space, interest_rate, n_t, n_s, sigma)


def option_initial_conditions(strike_price, put, x, n, s_star):
    """
    This function returns the initial conditions for the option
    :param strike_price: float The strike price of the option
    :param put: bool True if the option is a put option, False if it is a call option
    :param x: int The current x value
    :param n: int The total number of x values
    :param s_star: float The maximum value of the stock price
    :return: float The value of the option when the time is 0
    """

    # Todo: adjust this function to return the correct initial conditions for boundaries

    s = x * s_star / n
    if put:
        return max(strike_price - s, 0)
    else:
        return max(s - strike_price, 0)


def explicit_black_scholes(strike_price, put, boundary, time, s_star, interest_rate, n_t, n_s, sigma):
    h_t = time / n_t
    h_s = s_star / n_s

    lambda_0 = h_t / h_s ** 2

    if lambda_0 > 0.5:
        raise ValueError(f"lambda_0 must be less than 0.5, now it is {lambda_0}, h_t = {h_t}, h_x = {h_s}")

    # Add initial conditions to the grid
    solution = np.zeros((n_t, n_s))
    for i in range(n_s):
        solution[0, i] = option_initial_conditions(strike_price, put, i, n_s, s_star)

    # Add the boundary conditions to the grid
    solution[:, 0] = boundary[0]
    solution[:, -1] = boundary[1]

    for t in range(1, n_t):
        for i in range(1, n_s - 1):
            solution[t, i] = solution[t - 1, i - 1] * (h_t / 2.) * (sigma ** 2 * i ** 2 - interest_rate * i) + \
                             solution[t - 1, i] * (1 - h_t * (sigma ** 2) * (i ** 2) - h_t * interest_rate) + \
                             solution[t - 1, i + 1] * (h_t / 2.) * (sigma ** 2 * i ** 2 + interest_rate * i)

    return solution


def crank_nicolson_black_scholes(strike_price, put, boundary, time, s_star, interest_rate, n_t, n_s, sigma):
    h_t = time / n_t
    h_s = s_star / n_s

    # Add initial conditions to the grid
    solution = np.zeros((n_t, n_s))
    for i in range(n_s):
        solution[0, i] = option_initial_conditions(strike_price, put, i, n_s, s_star)

    # Add the boundary conditions to the grid
    solution[:, 0] = boundary[0]
    solution[:, -1] = boundary[1]

    left_matrix = np.zeros((n_s - 2, n_s - 2))
    right_matrix = np.zeros((n_s, n_s - 2))

    a_1 = None
    c_end = None

    for i in range(0, n_s - 2):
        i = i + 1
        a_i = - sigma ** 2 * h_t * i ** 2 / 4 + interest_rate * h_t * i / 4
        b_i = 1 + sigma ** 2 * h_t * i ** 2 / 2 + interest_rate * h_t / 2
        c_i = - sigma ** 2 * h_t * i ** 2 / 4 - interest_rate * h_t * i / 4
        d_i = 1 - sigma ** 2 * h_t * i ** 2 / 2 - interest_rate * h_t / 2


        i = i - 1

        if i != 0:
            left_matrix[i, i - 1] = a_i
        else:
            a_1 = a_i
        left_matrix[i, i] = b_i
        if i != n_s - 3:
            left_matrix[i, i + 1] = c_i
        else:
            c_end = c_i

        right_matrix[i:i + 3, i] = [-a_i, d_i, -c_i]

    right_matrix = right_matrix.transpose()



    # Apply the black-scholes equation
    for t in range(1, n_t):
        y = solution[t - 1, :]
        a = np.matmul(right_matrix, y)

        a[0] -= a_1 * solution[t, 0]
        a[-1] -= c_end * solution[t, -1]

        solution[t, 1:-1] = np.linalg.solve(left_matrix, a)

    assert np.all(np.isclose(solution[:, 0], boundary[0], atol=1e-3))
    assert np.all(np.isclose(solution[:, -1], boundary[1], atol=1e-3))

    return solution


if __name__ == '__main__':
    r = 0.06
    sigma = 0.3
    t_end = 1.
    strike = 10.
    s_star = 15.
    n_s = 50
    n_t = 50
    put = True

    boundary = [0, s_star - strike]
    if put:
        boundary = [strike, 0]

    result = apply_black_scholes(strike, put, boundary, t_end, s_star, r, n_t, n_s, sigma, crank_nicolson_black_scholes)

    # Plot initial conditions
    for i in range(0, n_t, 10):
        plt.plot(np.linspace(0, s_star, n_s), result[i, :])
        plt.show()

    # Plot the result as a plane
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    X = np.linspace(0, s_star, n_s)
    Y = np.linspace(0, t_end, n_t)
    X, Y = np.meshgrid(X, Y)
    #surf = ax.plot_surface(X, Y, result, cmap=cm.coolwarm,
    #                       linewidth=0, antialiased=False)

    #plt.show()

    # Plot the result as a gif of a rotating plane
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    def animate(i):
        ax.clear()
        ax.plot_surface(X, Y, result, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)
        ax.view_init(30, i)

    anim = animation.FuncAnimation(fig, animate, frames=360, interval=20, blit=False)
    anim.save('crank_nicolson_black_scholes.gif', writer='imagemagick', fps=30)
