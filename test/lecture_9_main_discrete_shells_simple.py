# -*- coding: utf-8 -*-
"""Lecture 9 - Main Discrete Shells Simple.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1el8_5wMuS4nNtu-hPZpMVVKzrmplcgn9

# Discrete Elastic Shells: Simple Example with Four Nodes
Copyright M. Khalid Jawed (khalidjm@seas.ucla.edu). License: CC BY-NC
You should use this code at your own risk.

#Load Libraries
"""

# from IPython.display import clear_output
import numpy as np
import matplotlib.pyplot as plt
import sys


"""#Miscellaneous Functions: signedAngle, rotateAxisAngle, parallel_transport, crossMat"""

def signedAngle(u = None,v = None,n = None):
    # This function calculates the signed angle between two vectors, "u" and "v",
    # using an optional axis vector "n" to determine the direction of the angle.
    #
    # Parameters:
    #   u: numpy array-like, shape (3,), the first vector.
    #   v: numpy array-like, shape (3,), the second vector.
    #   n: numpy array-like, shape (3,), the axis vector that defines the plane
    #      in which the angle is measured. It determines the sign of the angle.
    #
    # Returns:
    #   angle: float, the signed angle (in radians) from vector "u" to vector "v".
    #          The angle is positive if the rotation from "u" to "v" follows
    #          the right-hand rule with respect to the axis "n", and negative otherwise.
    #
    # The function works by:
    # 1. Computing the cross product "w" of "u" and "v" to find the vector orthogonal
    #    to both "u" and "v".
    # 2. Calculating the angle between "u" and "v" using the arctan2 function, which
    #    returns the angle based on the norm of "w" (magnitude of the cross product)
    #    and the dot product of "u" and "v".
    # 3. Using the dot product of "n" and "w" to determine the sign of the angle.
    #    If this dot product is negative, the angle is adjusted to be negative.
    #
    # Example:
    #   signedAngle(np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]))
    #   This would return a positive angle (π/2 radians), as the rotation
    #   from the x-axis to the y-axis is counterclockwise when viewed along the z-axis.
    w = np.cross(u,v)
    angle = np.arctan2( np.linalg.norm(w), np.dot(u,v) )
    if (np.dot(n,w) < 0):
        angle = - angle

    return angle

def mmt(matrix):
    return matrix + matrix.T

"""# Hinge angle, its gradient, and Hessian"""

#          x2
#          /\
#         /  \
#      e1/    \e3
#       /  t0  \
#      /        \
#     /    e0    \
#   x0------------x1
#     \          /
#      \   t1   /
#       \      /
#      e2\    /e4
#         \  /
#          \/
#          x3
#
#  Edge orientation: e0,e1,e2 point away from x0
#                       e3,e4 point away from x1

def getTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0

    n0 = np.cross(m_e0, m_e1)
    n1 = np.cross(m_e2, m_e0)

    # Calculate the signed angle using the provided function
    theta = signedAngle(n0, n1, m_e0)

    return theta

# In the original code, there are probaly TWO sign errors in the expressions for m_h3 and m_h4.
# [Original code: % https://github.com/shift09/plates-shells/blob/master/src/bending.cpp]
# I indicated those two corrections by writing the word "CORRECTION" next
# to them.

def gradTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0
    m_e3 = x2 - x1
    m_e4 = x3 - x1

    m_cosA1 = np.dot(m_e0, m_e1) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_cosA2 = np.dot(m_e0, m_e2) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_cosA3 = -np.dot(m_e0, m_e3) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_cosA4 = -np.dot(m_e0, m_e4) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_sinA1 = np.linalg.norm(np.cross(m_e0, m_e1)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_sinA2 = np.linalg.norm(np.cross(m_e0, m_e2)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_sinA3 = -np.linalg.norm(np.cross(m_e0, m_e3)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_sinA4 = -np.linalg.norm(np.cross(m_e0, m_e4)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_nn1 = np.cross(m_e0, m_e3)
    m_nn1 = m_nn1 / np.linalg.norm(m_nn1)
    m_nn2 = -np.cross(m_e0, m_e4)
    m_nn2 = m_nn2 / np.linalg.norm(m_nn2)

    m_h1 = np.linalg.norm(m_e0) * m_sinA1
    m_h2 = np.linalg.norm(m_e0) * m_sinA2
    m_h3 = -np.linalg.norm(m_e0) * m_sinA3  # CORRECTION
    m_h4 = -np.linalg.norm(m_e0) * m_sinA4  # CORRECTION
    m_h01 = np.linalg.norm(m_e1) * m_sinA1
    m_h02 = np.linalg.norm(m_e2) * m_sinA2

    # Initialize the gradient
    gradTheta = np.zeros(12)

    gradTheta[0:3] = m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4
    gradTheta[3:6] = m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2
    gradTheta[6:9] = -m_nn1 / m_h01
    gradTheta[9:12] = -m_nn2 / m_h02

    return gradTheta

def hessTheta(x0, x1 = None, x2 = None, x3 = None):

    if np.size(x0) == 12:  # Allow another type of input where x0 contains all the info
      x1 = x0[3:6]
      x2 = x0[6:9]
      x3 = x0[9:12]
      x0 = x0[0:3]

    m_e0 = x1 - x0
    m_e1 = x2 - x0
    m_e2 = x3 - x0
    m_e3 = x2 - x1
    m_e4 = x3 - x1

    m_cosA1 = np.dot(m_e0, m_e1) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_cosA2 = np.dot(m_e0, m_e2) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_cosA3 = -np.dot(m_e0, m_e3) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_cosA4 = -np.dot(m_e0, m_e4) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_sinA1 = np.linalg.norm(np.cross(m_e0, m_e1)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e1))
    m_sinA2 = np.linalg.norm(np.cross(m_e0, m_e2)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e2))
    m_sinA3 = -np.linalg.norm(np.cross(m_e0, m_e3)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e3))
    m_sinA4 = -np.linalg.norm(np.cross(m_e0, m_e4)) / (np.linalg.norm(m_e0) * np.linalg.norm(m_e4))

    m_nn1 = np.cross(m_e0, m_e3)
    m_nn1 /= np.linalg.norm(m_nn1)
    m_nn2 = -np.cross(m_e0, m_e4)
    m_nn2 /= np.linalg.norm(m_nn2)

    m_h1 = np.linalg.norm(m_e0) * m_sinA1
    m_h2 = np.linalg.norm(m_e0) * m_sinA2
    m_h3 = -np.linalg.norm(m_e0) * m_sinA3
    m_h4 = -np.linalg.norm(m_e0) * m_sinA4
    m_h01 = np.linalg.norm(m_e1) * m_sinA1
    m_h02 = np.linalg.norm(m_e2) * m_sinA2

    # Gradient of Theta (as an intermediate step)
    grad_theta = np.zeros((12, 1))
    grad_theta[0:3] = (m_cosA3 * m_nn1 / m_h3 + m_cosA4 * m_nn2 / m_h4).reshape(-1, 1)
    grad_theta[3:6] = (m_cosA1 * m_nn1 / m_h1 + m_cosA2 * m_nn2 / m_h2).reshape(-1, 1)
    grad_theta[6:9] = (-m_nn1 / m_h01).reshape(-1, 1)
    grad_theta[9:12] = (-m_nn2 / m_h02).reshape(-1, 1)

    # Intermediate matrices for Hessian
    m_m1 = np.cross(m_nn1, m_e1) / np.linalg.norm(m_e1)
    m_m2 = -np.cross(m_nn2, m_e2) / np.linalg.norm(m_e2)
    m_m3 = -np.cross(m_nn1, m_e3) / np.linalg.norm(m_e3)
    m_m4 = np.cross(m_nn2, m_e4) / np.linalg.norm(m_e4)
    m_m01 = -np.cross(m_nn1, m_e0) / np.linalg.norm(m_e0)
    m_m02 = np.cross(m_nn2, m_e0) / np.linalg.norm(m_e0)

    # Hessian matrix components
    M331 = m_cosA3 / (m_h3 ** 2) * np.outer(m_m3, m_nn1)
    M311 = m_cosA3 / (m_h3 * m_h1) * np.outer(m_m1, m_nn1)
    M131 = m_cosA1 / (m_h1 * m_h3) * np.outer(m_m3, m_nn1)
    M3011 = m_cosA3 / (m_h3 * m_h01) * np.outer(m_m01, m_nn1)
    M111 = m_cosA1 / (m_h1 ** 2) * np.outer(m_m1, m_nn1)
    M1011 = m_cosA1 / (m_h1 * m_h01) * np.outer(m_m01, m_nn1)

    M442 = m_cosA4 / (m_h4 ** 2) * np.outer(m_m4, m_nn2)
    M422 = m_cosA4 / (m_h4 * m_h2) * np.outer(m_m2, m_nn2)
    M242 = m_cosA2 / (m_h2 * m_h4) * np.outer(m_m4, m_nn2)
    M4022 = m_cosA4 / (m_h4 * m_h02) * np.outer(m_m02, m_nn2)
    M222 = m_cosA2 / (m_h2 ** 2) * np.outer(m_m2, m_nn2)
    M2022 = m_cosA2 / (m_h2 * m_h02) * np.outer(m_m02, m_nn2)

    B1 = 1 / np.linalg.norm(m_e0) ** 2 * np.outer(m_nn1, m_m01)
    B2 = 1 / np.linalg.norm(m_e0) ** 2 * np.outer(m_nn2, m_m02)

    N13 = 1 / (m_h01 * m_h3) * np.outer(m_nn1, m_m3)
    N24 = 1 / (m_h02 * m_h4) * np.outer(m_nn2, m_m4)
    N11 = 1 / (m_h01 * m_h1) * np.outer(m_nn1, m_m1)
    N22 = 1 / (m_h02 * m_h2) * np.outer(m_nn2, m_m2)
    N101 = 1 / (m_h01 ** 2) * np.outer(m_nn1, m_m01)
    N202 = 1 / (m_h02 ** 2) * np.outer(m_nn2, m_m02)

    # Initialize Hessian of Theta
    hess_theta = np.zeros((12, 12))

    hess_theta[0:3, 0:3] = mmt(M331) - B1 + mmt(M442) - B2
    hess_theta[0:3, 3:6] = M311 + M131.T + B1 + M422 + M242.T + B2
    hess_theta[0:3, 6:9] = M3011 - N13
    hess_theta[0:3, 9:12] = M4022 - N24
    hess_theta[3:6, 3:6] = mmt(M111) - B1 + mmt(M222) - B2
    hess_theta[3:6, 6:9] = M1011 - N11
    hess_theta[3:6, 9:12] = M2022 - N22
    hess_theta[6:9, 6:9] = -mmt(N101)
    hess_theta[9:12, 9:12] = -mmt(N202)

    # Make the Hessian symmetric
    hess_theta[3:6, 0:3] = hess_theta[0:3, 3:6].T
    hess_theta[6:9, 0:3] = hess_theta[0:3, 6:9].T
    hess_theta[9:12, 0:3] = hess_theta[0:3, 9:12].T
    hess_theta[6:9, 3:6] = hess_theta[3:6, 6:9].T
    hess_theta[9:12, 3:6] = hess_theta[3:6, 9:12].T

    return hess_theta

"""# Stretching energy for a shell, it's gradient, and Hessian"""

def gradEs_hessEs(node0 = None,node1 = None,l_k = None,EA = None):

# Inputs:
# node0: 1x3 vector - position of the first node
# node1: 1x3 vector - position of the last node

# l_k: reference length (undeformed) of the edge
# EA: scalar - stretching stiffness - Young's modulus times area

# Outputs:
# dF: 6x1  vector - gradient of the stretching energy between node0 and node 1.
# dJ: 6x6 vector - hessian of the stretching energy between node0 and node 1.

    ## Gradient of Es
    edge = node1 - node0

    edgeLen = np.linalg.norm(edge)
    tangent = edge / edgeLen
    epsX = edgeLen / l_k - 1
    dF_unit = EA * tangent * epsX
    dF = np.zeros((6))
    dF[0:3] = - dF_unit
    dF[3:6] = dF_unit

    ## Hessian of Es
    Id3 = np.eye(3)
    M = EA * ((1 / l_k - 1 / edgeLen) * Id3 + 1 / edgeLen * ( np.outer( edge, edge ) ) / edgeLen ** 2)

    dJ = np.zeros((6,6))
    dJ[0:3,0:3] = M
    dJ[3:6,3:6] = M
    dJ[0:3,3:6] = - M
    dJ[3:6,0:3] = - M
    return dF,dJ

"""# Bending energy for a shell, it's gradient, and Hessian"""

def getEb_Shell(x0, x1=None, x2=None, x3=None, theta_bar=0, kb=1.0):
    """
    Compute the bending energy for a shell.

    Returns:
    E (scalar): Bending energy.
    """
    # Allow another type of input where x0 contains all the information
    if np.size(x0) == 12:
        x1 = x0[3:6]
        x2 = x0[6:9]
        x3 = x0[9:12]
        x0 = x0[:3]

    # Compute theta, gradient, and Hessian
    theta = getTheta(x0, x1, x2, x3)  # Replace with your getTheta function in Python
    grad = gradTheta(x0, x1, x2, x3)  # Replace with your gradTheta function in Python

    # E = 0.5 * kb * (theta-thetaBar)^2
    E = 0.5 * kb * (theta - theta_bar) ** 2

    return E

def gradEb_hessEb_Shell(x0, x1=None, x2=None, x3=None, theta_bar=0, kb=1.0):
    """
    Compute the gradient and Hessian of the bending energy for a shell.

    Parameters:
    x0 (array): Can either be a 3-element array (single point) or a 12-element array.
    x1, x2, x3 (arrays): Optional, 3-element arrays specifying points.
    theta_bar (float): Reference angle.
    kb (float): Bending stiffness.

    Returns:
    dF (array): Gradient of the bending energy.
    dJ (array): Hessian of the bending energy.
    """
    # Allow another type of input where x0 contains all the information
    if np.size(x0) == 12:
        x1 = x0[3:6]
        x2 = x0[6:9]
        x3 = x0[9:12]
        x0 = x0[:3]

    # Compute theta, gradient, and Hessian
    theta = getTheta(x0, x1, x2, x3)  # Replace with your getTheta function in Python
    grad = gradTheta(x0, x1, x2, x3)  # Replace with your gradTheta function in Python

    # E = 0.5 * kb * (theta-thetaBar)^2
    # F = dE/dx = 2 * (theta-thetaBar) * gradTheta
    dF = 0.5 * kb * (2 * (theta - theta_bar) * grad)

    # E = 0.5 * kb * (theta-thetaBar)^2
    # F = 0.5 * kb * (2 (theta-thetaBar) d theta/dx)
    # J = dF/dx = 0.5 * kb * [ 2 (d theta / dx) transpose(d theta/dx) +
    #       2 (theta-thetaBar) (d^2 theta/ dx^2 ) ]
    hess = hessTheta(x0, x1, x2, x3)  # Replace with your hessTheta function in Python
    dJ = 0.5 * kb * (2 * np.outer(grad, grad) + 2 * (theta - theta_bar) * hess)

    return dF, dJ

"""#Plot the shell"""

# Function to set equal aspect ratio for 3D plots
def set_axes_equal(ax):
    """
    Set equal aspect ratio for a 3D plot in Matplotlib.
    This function adjusts the limits of the plot to make sure
    that the scale is equal along all three axes.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max(x_range, y_range, z_range)

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2])

def plotShell(x0, ctime):

  x1 = x0[3:6]
  x2 = x0[6:9]
  x3 = x0[9:12]
  x0 = x0[0:3]

  fig = plt.figure(1)
  # clear_output()
  plt.clf()  # Clear the figure
  ax = fig.add_subplot(111, projection='3d')

  # Plot nodes
  X = np.array([x0[0], x1[0], x2[0], x0[0], x3[0], x1[0]])
  Y = np.array([x0[1], x1[1], x2[1], x0[1], x3[1], x1[1]])
  Z = np.array([x0[2], x1[2], x2[2], x0[2], x3[2], x1[2]])
  ax.plot3D(X, Y, Z, 'ko-')

  # Plot the first node with a red triangle
  ax.plot3D([X[0]], [Y[0]], [Z[0]], 'r^')

  # Set the title with current time
  ax.set_title(f't={ctime:.2f}')

  # Set axes labels
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')

  # Set equal scaling and a 3D view
  set_axes_equal(ax)
  plt.draw()  # Force a redraw of the figure

  plt.show()

"""# Objective function"""

def objfun(qGuess, q0, u, freeIndex, dt, tol,
                massVector, mMat,
                ks, lk, edges,
                kb, thetaBar, hinges,
                Fg):

  q = qGuess
  ndof = len(q)
  iter = 0
  error = 10 * tol

  while error > tol:

    # Calculating Bending
    Fb = np.zeros(ndof)
    Jb = np.zeros((ndof,ndof))
    for kHinge in range(hinges.shape[0]):
      node0 = hinges[kHinge,0]
      node1 = hinges[kHinge,1]
      node2 = hinges[kHinge,2]
      node3 = hinges[kHinge,3]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      x2 = q[3*node2:3*node2+3]
      x3 = q[3*node3:3*node3+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2,
             3*node2, 3*node2 + 1, 3*node2 + 2,
             3*node3, 3*node3 + 1, 3*node3 + 2]
      dF, dJ = gradEb_hessEb_Shell(x0, x1, x2, x3, thetaBar[kHinge], kb)
      Fb[ind] = Fb[ind] - dF
      Jb[np.ix_(ind,ind)] = Jb[np.ix_(ind,ind)] - dJ

    # Calculating stretching
    Fs = np.zeros(ndof)
    Js = np.zeros((ndof,ndof))
    for kEdge in range(edges.shape[0]):
      node0 = edges[kEdge,0]
      node1 = edges[kEdge,1]
      x0 = q[3*node0:3*node0+3]
      x1 = q[3*node1:3*node1+3]
      ind = [3*node0, 3*node0 + 1, 3*node0 + 2,
             3*node1, 3*node1 + 1, 3*node1 + 2]
      dF, dJ = gradEs_hessEs(x0, x1, lk[kEdge], ks[kEdge])
      Fs[ind] = Fs[ind] - dF
      Js[np.ix_(ind,ind)] = Js[np.ix_(ind,ind)] - dJ

    # Calculating total force
    F = Fg + Fb + Fs # Viscous forces can sometimes be useful
    JForces = Js + Jb

    JForces_free = JForces[np.ix_(freeIndex, freeIndex)]
    mMat_free = mMat[np.ix_(freeIndex, freeIndex)]
    f = massVector/dt * ( (q-q0)/dt - u ) - F
    J = mMat / dt ** 2 - JForces

    # print(JForces_free)
    # print(mMat_free)

    f_free = f[freeIndex]
    J_free = J[np.ix_(freeIndex, freeIndex)]
    dq_free = np.linalg.solve(J_free, f_free)

    # print(J_free)

    q[freeIndex] = q[freeIndex] - dq_free

    error = np.sum( np.abs(f_free) )
    iter += 1

    # print('Iter = %d' % iter)
    # print('Error = %f' % error)
  print(q[freeIndex])
  u = (q - q0) / dt

  return q, u

"""#Main Discrete Shells

**Degrees of freedom and nodes**

We are dividing an elastic shell  into nv nodes, which corresponds to a DOF vector of size $$ndof=3nv.$$.
"""

# Inputs


# x1 = np.array([0.01, 0.0, 0.0])
# x2 = np.array([0.005, 0.01, 0.0])
# x3 = np.array([0.005, -0.01, 0.0])

# edges = np.array([(0,1), (1,2), (0,2), (0,3), (1,3)])
# hinges = np.array([(0,1,2,3)])


# x0 = np.array([0.0, 0.0, 0.0]) 
# x1 = np.array([0.1, 0.0, 0.0])
# x2 = np.array([0.05, 0.0, -0.0866025388])
# x3 = np.array([0.15, 0.0, -0.0866025388])


# q0 = np.array([0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.05, 0.0, -0.0866025388, 0.15, 0.0, -0.0866025388])
# edges = np.array([(2, 1), (1, 0), (0, 2), (1, 3), (2, 3)])
# hinges = np.array([(2, 1, 0, 3)])
# nv = 4 # Number of vertices
np.set_printoptions(formatter={'float': '{: .8e}'.format})



q0 = np.array([ 0.00000000,  0.00000000,  0.00000000,  0.10000000,  0.00000000,  0.00000000,  0.20000000,  0.00000000,  0.00000000,  0.05000000,  0.00000000, -0.08660254,  0.15000000,  0.00000000, -0.08660254,  0.25000000,  0.00000000, -0.08660254,  0.00000000,  0.00000000, -0.17320508,  0.10000000,  0.00000000, -0.17320508,  0.20000000,  0.00000000, -0.17320508])
edges = np.array([(3, 1), (1, 0), (0, 3), (4, 2), (2, 1), (1, 4), (2, 5), (3, 6), (7, 4), (4, 3), (3, 7), (8, 5), (5, 4), (4, 8), (6, 7), (7, 8)])
hinges = np.array([(3, 1, 0, 4), (4, 2, 1, 5), (1, 4, 2, 3), (3, 7, 4, 6), (3, 4, 1, 7), (4, 8, 5, 7), (4, 7, 8, 3), (4, 5, 2, 8)])
nv = q0.shape[0] // 3 # number of vertices

ndof = 3 * nv # number of DOFs

fixedIndex = np.arange(0,9)
freeIndex = np.arange(9,ndof)

orig_stdout = sys.stdout

nv_side = int(np.sqrt(nv)) - 1
filename = f"./test/py{nv_side}.out.txt"
f = open(filename, 'w')
sys.stdout = f

"""**Create Edges and Hinges**"""





"""**Elastic stiffness parameters**
We will compute the elastic stiffness values
"""

Y = 1.0e7 # Pascals
h = 0.001 # meter (thickness)

# Bending stiffness
kb = 2.0/np.sqrt(3.0) * Y * h**3.0 / 12.0

# Stretching stiffness and reference/undeformed length (lk)
lk = np.zeros(edges.shape[0]) # edges.shape[0] is number of edges
ks = np.zeros(edges.shape[0])
for kEdge in range(edges.shape[0]):
  node0 = edges[kEdge,0] # node number
  node1 = edges[kEdge,1] # node numer
  n0 = q0[3*node0:3*node0+3]
  n1 = q0[3*node1:3*node1+3]
  lk[kEdge] = np.linalg.norm(n1-n0)
  ks[kEdge] = np.sqrt(3.0)/2.0 * Y * h * (lk[kEdge]) ** 2.0

"""**Time parameters**"""

totalTime = 5 # Simulation time
dt = 0.01 # Time step size

# Tolerance
tol = kb / (0.01) * 1.0e-3

"""**Mass vector and matrix**"""

totalM = 0.01 # 10 g is total mass
dm = totalM / nv # mass per node: simplification but OK here
massVector = np.zeros(ndof)
for kNode in range(nv):
  massVector[3*kNode:3*kNode+3] = dm * np.ones(3)
mMat = np.diag(massVector)

"""**External Force**
Since this does not change with time, we define it here instead of within the time loop.
"""

g = np.array([0.0, -9.81, 0.0])
Fg = np.zeros(ndof)
for kNode in range(nv):
  Fg[3*kNode:3*kNode+3] = dm * g

"""**Initial DOF vector**"""

# Initial position: q0 is already define
u = np.zeros(ndof) # initial velocity

"""**Compute natural curvature (theta_bar)**"""

# ss = np.concatenate((x2, x1, x0, x3))
# thetaBar = getTheta(ss) # Just 0 for plates
thetaBar = np.zeros(hinges.shape[0])
for kHinge in range(hinges.shape[0]):
  node0 = hinges[kHinge,0]
  node1 = hinges[kHinge,1]
  node2 = hinges[kHinge,2]
  node3 = hinges[kHinge,3]
  x0 = q0[3*node0:3*node0+3]
  x1 = q0[3*node1:3*node1+3]
  x2 = q0[3*node2:3*node2+3]
  x3 = q0[3*node3:3*node3+3]
  thetaBar[kHinge] = getTheta(x0, x1, x2, x3)

"""**Set up boundary conditions**
We define the fixed and free DOFs. If the simulation asks that the boundary conditions vary with time (e.g., someone is holding a rod and all on a sudden drops it), we have to define it later within the time stepping loop.
"""

# freeIndex = np.arange(9,ndof)

"""**Time stepping loop**"""

Nsteps = round(totalTime / dt)
ctime = 0.0 # Current time
endZ = np.zeros(Nsteps) # Store z-coordinate of last node i.e. tip deflection

for timeStep in range(Nsteps):
  print('Time: %f' % ctime)

  # Call objective function
  qGuess = q0.copy()
  q, u = objfun(qGuess, q0, u, freeIndex, dt, tol,
                massVector, mMat,
                ks, lk, edges,
                kb, thetaBar, hinges,
                Fg)
  ctime += dt # Increment time

  q0 = q.copy()
  endZ[timeStep] = q[-1] # Store z-coordinate of last node

  # if timeStep % 100 == 0:
  #   plotShell(q, ctime)

sys.stdout = orig_stdout
f.close()

# Visualize
plt.figure(2)
time_array = np.arange(1, Nsteps+1) * dt
plt.plot(time_array, endZ, 'ro-')
plt.xlabel('TIme, t [sec]')
plt.ylabel('z-coord of last node, $\\delta_z$ [m]')
plt.title('Tip deflection vs. time')
plt.show()