"""
This is a script to test different versions of the finite difference method to solve the heat equation, 
and later the Black-Scholes equation
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def apply_heat_equation(initial, boundary, time, space, time_step, space_step, heat_equation):
    return heat_equation(initial, boundary, time, space, time_step, space_step)


def explicit_heat_equation(initial, boundary, time, space, time_step, space_step):
    """
    This is the explicit method, which is the fastest, but also the least accurate (according to Github Copilot). It has
    a stability limit of lambda ≤ 0.5 and an error of O(h_t + O(h_x^2)
    :param initial: The initial conditions, as a function of x and n
    :param boundary: The boundary conditions, as a list of two values
    :param time: The time interval
    :param space: The space interval
    :param time_step: The time step size as a float
    :param space_step: The space step size as a float
    :return: A matrix of the solution
    """


    lambda_0 = time_step / space_step ** 2
    lambda_1 = 1 - 2 * lambda_0

    n = int(space / space_step)
    m = int(time / time_step)

    if lambda_0 > 0.5:
        raise ValueError(f"lambda_0 must be less than 0.5, now it is {lambda_0}, h_t = {time_step}, h_x = {space_step}")

    # Add initial conditions to the grid
    solution = np.zeros((m, n))
    for i in range(n):
        solution[0, i] = initial(i, n)

    # Add the boundary conditions to the grid
    solution[:, 0] = boundary[0]
    solution[:, -1] = boundary[1]

    # Apply the heat equation
    for t in range(1, m):
        solution[t, 1:-1] = lambda_0 * solution[t - 1, 2:] + lambda_1 * solution[t - 1, 1:-1] + lambda_0 * solution[
                                                                                                           t - 1, :-2]

    return solution


def implicit_heat_equation(initial, boundary, time, space, time_step, space_step):
    """
    This is the implicit method, which is the most accurate, but also the slowest (according to Github Copilot). It has
    an error of O(h_x^2) + O(h_t) and is stable for any lambda_0
    :param initial: The initial conditions, as a function of x and n
    :param boundary: The boundary conditions, as a list of two values
    :param time: The time interval
    :param space: The space interval
    :param time_step: The time step size as a float
    :param space_step: The space step size as a float
    :return: A matrix of the solution
    """
    lambda_0 = -time_step / space_step ** 2
    lambda_1 = 1 - 2 * lambda_0

    n = int(space / space_step)
    m = int(time / time_step)

    # Add initial conditions to the grid
    solution = np.zeros((m, n))
    for i in range(n):
        solution[0, i] = initial(i, n-1)

    # Add the boundary conditions to the grid
    solution[:, 0] = boundary[0]
    solution[:, -1] = boundary[1]

    # Create 3 diagonals of the matrix
    diagonal_0 = np.ones(n - 3) * lambda_0
    diagonal_1 = np.ones(n - 2) * lambda_1
    diagonal_2 = np.ones(n - 3) * lambda_0

    # Combine the diagonals into a matrix
    matrix = np.diag(diagonal_0, -1) + np.diag(diagonal_1, 0) + np.diag(diagonal_2, 1)

    assert matrix.shape == (n - 2, n - 2)

    # Apply the heat equation
    for t in range(1, int(time / time_step)):
        y = solution[t - 1, 1:-1]
        y[0] += lambda_0 * boundary[0]
        y[-1] += lambda_0 * boundary[1]
        solution[t, 1:-1] = np.linalg.solve(matrix, y)

    assert solution.shape == (m, n)
    assert np.max(solution[:, 0]) == boundary[0]
    assert np.max(solution[:, -1]) == boundary[1]

    return solution


def crank_nicolson_heat_equation(initial, boundary, time, space, time_step, space_step):
    """
    This is the Crank-Nicolson method, which is a combination of the explicit and implicit method

    :param initial: The initial conditions, as a function of x and n
    :param boundary: The boundary conditions, as a list of two values
    :param time: The time interval
    :param space: The space interval
    :param time_step: The time step size as a float
    :param space_step: The space step size as a float
    :return: A matrix of the solution
    """
    lambda_0 = time_step / space_step ** 2

    n = int(space / space_step)
    m = int(time / time_step)

    # Add initial conditions to the grid
    solution = np.zeros((m, n))
    for i in range(n):
        solution[0, i] = initial(i, n-1)

    solution[:, 0] = boundary[0]
    solution[:, -1] = boundary[1]

    # Create 3 diagonals of the a_matrix
    diagonal_0 = np.ones(n - 2) * (lambda_0 / 2.)
    diagonal_1 = np.ones(n - 2) * (1 - lambda_0)
    diagonal_2 = np.ones(n - 2) * (lambda_0 / 2.)

    # Combine the diagonals into a matrix
    matrix_a = np.zeros((n - 2, n))

    matrix_a[:, 0:-2] += np.diag(diagonal_0, 0)  # np.diag(diagonal_1, 1) + np.diag(diagonal_2, 0)
    matrix_a[:, 2:] += np.diag(diagonal_0, 0)
    matrix_a[:, 1:-1] += np.diag(diagonal_1, 0)

    assert matrix_a.shape == (n - 2, n)

    # Create 3 diagonals of the b_matrix
    diagonal_0 = np.ones(n - 3) * (-lambda_0 / 2.)
    diagonal_1 = np.ones(n - 2) * (1 + lambda_0)
    diagonal_2 = np.ones(n - 3) * (-lambda_0 / 2.)

    # Combine the diagonals into a matrix
    matrix_b = np.diag(diagonal_0, -1) + np.diag(diagonal_1, 0) + np.diag(diagonal_2, 1)

    assert matrix_b.shape == (n - 2, n - 2)

    # Apply the heat equation
    for t in range(1, int(time / time_step)):
        y = solution[t - 1, :]
        a = np.matmul(matrix_a, y)

        a[0] += lambda_0 * solution[t, 0]
        a[-1] += lambda_0 * solution[t, -1]

        solution[t, 1:-1] = np.linalg.solve(matrix_b, a)

    assert solution.shape == (m, n)
    assert np.max(solution[:, 0]) == boundary[0]
    assert np.max(solution[:, -1]) == boundary[1]

    return solution

def calculate_error(function, time, space, time_step, space_step, result):

    n = int(space / space_step)
    m = int(time / time_step)

    real_solution = np.zeros((m, n))

    for i in range(n):
        for j in range(m):
            real_solution[j, i] = function(i, n-1, j, m-1)

    assert np.max(real_solution[:, 0]) == 0.
    assert np.max(real_solution[:, -1]) == 0.

    error = np.abs(result - real_solution)

    return error


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # Define the initial conditions
    def initial_function(x_i, x_max, t_j=0., t_max=1.):
        if x_i in [0, x_max]:
            return 0.0

        x = x_i / x_max
        t = t_j / t_max
        return np.sin(np.pi * x) * np.exp(t * -np.pi ** 2)


    # Define the boundary conditions
    boundary_conditions = [0., 0.]

    time = 1
    space = 1

    time_step = 0.01
    space_step = 0.0001

    method_dict = {
        "Crank Nicolson": crank_nicolson_heat_equation,
        "Implicit Method": implicit_heat_equation,
        #"Explicit": explicit_heat_equation
    }

    method_error_dict = {
        "Crank Nicolson": [],
        "Implicit Method": [],
        #"Explicit": []
    }
    time_steps = []

    for multiplier in range(10):

        time_step = (1./128.) * (multiplier + 1)

        time_steps.append(time_step)

        for method in method_dict:
            result = apply_heat_equation(initial_function,
                                         boundary_conditions,
                                         time,
                                         space,
                                         time_step,
                                         space_step,
                                         method_dict[method])
            error = calculate_error(initial_function, time, space, time_step, space_step, result)

            error_list = method_error_dict[method]

            error_list.append(np.max(error))

            # for i in range(0, int(time / time_step), 5):
            #     plt.plot(np.linspace(0, space, int(space / space_step)), error[i, :])
            #
            # plt.title(f"Errors with {method}, h_t = {time_step}, h_x = {space_step}")
            # plt.show()
            #
            # # Plot the result as a plane
            # fig = plt.figure()
            # ax = fig.gca(projection='3d')
            # X = np.linspace(0, space, int(space / space_step))
            # Y = np.linspace(0, time, int(time / time_step))
            # X, Y = np.meshgrid(X, Y)
            # surf = ax.plot_surface(X, Y, error, cmap=cm.coolwarm,
            #                        linewidth=0, antialiased=False)
            #
            # plt.title(f"Errors with {method}, h_t = {time_step}, h_x = {space_step}")
            #
            # plt.show()


    for method in method_error_dict:
        plt.cla()
        plt.scatter(time_steps, method_error_dict[method])
        plt.title(f"Max error with h_x = {space_step} with method {method}")
        plt.xlabel("h_t")
        plt.ylabel("maximum error")
        plt.show()


    result = apply_heat_equation(initial_function,
                                 boundary_conditions,
                                 time,
                                 space,
                                 time_step,
                                 space_step,
                                 #implicit_heat_equation)
                                 crank_nicolson_heat_equation)



    # # Plot initial conditions
    # for i in range(0, int(time / time_step), 500):
    #     plt.plot(np.linspace(0, space, int(space / space_step)), result[i, :])
    # plt.show()
    #
    # # Plot the result as a plane
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # X = np.linspace(0, space, int(space / space_step))
    # Y = np.linspace(0, time, int(time / time_step))
    # X, Y = np.meshgrid(X, Y)
    # surf = ax.plot_surface(X, Y, result, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    #
    # plt.show()

    # Calculate the error
    error = calculate_error(initial_function, time, space, time_step, space_step, result)

    # Plot initial conditions
    plt.cla()

    for i in range(0, int(time / time_step), 5):
        plt.plot(np.linspace(0, space, int(space / space_step)), error[i, :])
    plt.show()

    # Plot the result as a plane
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0, space, int(space / space_step))
    Y = np.linspace(0, time, int(time / time_step))
    X, Y = np.meshgrid(X, Y)
    surf = ax.plot_surface(X, Y, error, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    plt.show()
