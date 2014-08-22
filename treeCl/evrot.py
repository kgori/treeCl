from __future__ import print_function
# third party
import numpy as np

# treeCl
from evrot_extensions import build_Uab, sum_dJ


def buildA(X, U1, Vk, U2):
    # A = X*U1*Vk*U2
    A = X.dot(U1.dot(Vk.dot(U2)))

    return A


def gradU(theta, k, ik, jk, dim):
    V = np.zeros((dim, dim), dtype=np.float)
    i = ik[k]
    j = jk[k]
    tt = theta[k]
    sin = np.sin(tt)
    cos = np.cos(tt)

    V[[i, i, j, j], [i, j, i, j]] = [-sin, cos, -cos, -sin]

    return V


def rotate_givens(X, theta, ik, jk, angle_num, dim):
    G = build_Uab(theta, 0, angle_num - 1, ik, jk, dim)
    Y = X.dot(G)
    return Y


def evqualitygrad(X, theta, ik, jk, angle_num, angle_index, dim, ndata):
    V = gradU(theta, angle_index, ik, jk, dim)

    U1 = build_Uab(theta, 0, angle_index - 1, ik, jk, dim)
    U2 = build_Uab(theta, angle_index + 1, angle_num - 1, ik, jk, dim)

    A = buildA(X, U1, V, U2)

    Y = rotate_givens(X, theta, ik, jk, angle_num, dim)

    # Find max of each row
    r_ndata = range(ndata)
    Y_sq = Y ** 2
    max_index = Y_sq.argmax(axis=1)
    max_values = Y[r_ndata, max_index]

    max_A_values = A[r_ndata, max_index]
    mv_sq = max_values ** 2
    mv_cb = max_values * mv_sq
    A_x_Y = A * Y

    # Compute gradient
    dJ = sum_dJ(A_x_Y, Y_sq, mv_sq, mv_cb, max_A_values, dim, ndata)

    return dJ


def evqual(X, ndata, dim):
    Xsquare = X * X
    max_values = Xsquare.max(axis=1)
    # max_index = Xsquare.argmax(axis=1)

    # Compute cost
    Xsq_div_maxvals = Xsquare / max_values.reshape((ndata, 1))
    sums = np.sum(Xsq_div_maxvals)
    J = 1 - (sums / ndata - 1) / dim

    return float(J)


def cluster_assign(X):  # , ik, jk, dim, ndata):

    (ndata, dim) = X.shape
    Xsq = X * X
    # max_values = np.max(Xsq, axis=1)
    max_index = np.argmax(Xsq, axis=1)
    cluster_count = np.bincount(max_index, minlength=dim)
    cluster_cell_array = [np.array([0] * count) for count in cluster_count]

    for j in xrange(dim):
        cluster = cluster_cell_array[j]
        cind = 0
        for i in xrange(ndata):
            if max_index[i] == j:
                cluster[cind] = i + 1
                cind += 1
    return cluster_cell_array


def test(X):
    (ndata, dim) = X.shape
    (ik, jk) = np.triu_indices(dim, k=1)
    angle_num = len(ik)

    theta = np.arange(0.0, angle_num / 10., 0.1)
    Q = evqual(X, ndata, dim)
    print('Q = {0}'.format(Q))
    dQ = evqualitygrad(X, theta, ik, jk, 45, 5, 10, 40)
    print('x:', X)
    print('theta:', theta)
    print('ik and jk', ik, jk)
    print('Q:', Q)
    print('dQ:', dQ)


def main(X, max_iter=200):
    (ndata, dim) = X.shape
    (ik, jk) = np.triu_indices(dim, k=1)
    angle_num = len(ik)

    theta = np.array([0.0] * angle_num)
    theta_new = np.array([0.0] * angle_num)

    Q = evqual(X, ndata, dim)
    Q_old1 = Q
    Q_old2 = Q
    alpha = 1
    for iteration in range(max_iter):

        for d in range(angle_num):
            dQ = evqualitygrad(X, theta, ik, jk, angle_num, d, dim, ndata)
            theta_new[d] = theta[d] - alpha * dQ
            Xrot = rotate_givens(X, theta_new, ik, jk, angle_num, dim)
            Q_new = evqual(Xrot, ndata, dim)
            if Q_new > Q:
                theta[d] = theta_new[d]
                Q = Q_new
            else:
                theta_new[d] = theta[d]


        # Stopping criterion
        if iteration > 1:
            if Q - Q_old2 < 0.001:
                break

        Q_old2 = Q_old1
        Q_old1 = Q

    Xrot = rotate_givens(X, theta_new, ik, jk, angle_num, dim)
    clusts = cluster_assign(Xrot)

    return (clusts, Q, Xrot)


if __name__ == "__main__":
    mat = np.array([[-0.1246, 0.0751, -0.1526, 0.1311, -0.0268, 0.2671, -0.2424, 0.0099, -0.0636, 0.4],
                    [-0.1649, 0.1295, 0.2068, -0.1068, 0.0568, 0.1167, 0.0181, -0.3972, 0.1179, 0.0022],
                    [-0.1601, 0.1252, 0.1979, -0.0785, 0.0398, 0.0794, -0.0002, -0.4268, 0.1202, 0.0045],
                    [-0.1045, -0.0736, -0.0201, 0.2591, 0.1855, 0.0789, 0.1621, -0.0451, -0.099, 0.2491],
                    [-0.1179, -0.1618, 0.0049, 0.1046, 0.3152, -0.0138, 0.0365, -0.1165, -0.4487, -0.088],
                    [-0.1776, 0.1405, 0.2136, 0.0453, -0.0422, -0.2152, -0.2001, 0.078, -0.0746, -0.0087],
                    [-0.1745, -0.2466, 0.01, 0.0413, 0.2065, -0.0396, -0.0783, 0.0477, 0.1563, -0.068],
                    [-0.1726, 0.114, -0.2056, 0.143, -0.0696, 0.1301, 0.0463, 0.0072, 0.047, -0.3139],
                    [-0.1771, -0.2524, 0.0113, -0.0141, 0.057, -0.0327, -0.082, 0.0601, 0.1981, -0.0917],
                    [-0.1401, 0.9741, -0.1979, -0.2004, 0.0791, -0.2961, 0.2447, -0.0143, 0.0115, 0.28],
                    [-0.1579, -0.2183, 0.0099, -0.0256, -0.0764, -0.0023, -0.0215, 0.0748, 0.3047, 0.1203],
                    [-0.1439, -0.1965, 0.0055, 0.1359, 0.3949, -0.0216, 0.0187, -0.0852, -0.3276, -0.0504],
                    [-0.1514, 0.1016, -0.2046, 0.0418, -0.0264, 0.2631, -0.3644, 0.0204, -0.0674, 0.2416],
                    [-0.1702, 0.1176, -0.2394, -0.125, 0.0332, -0.0818, 0.0186, -0.0051, 0.0115, -0.2969],
                    [-0.1665, 0.1294, 0.1837, 0.1455, -0.0931, -0.2642, -0.1562, -0.0791, -0.019, 0.0026],
                    [-0.1497, 0.086, -0.0804, 0.409, -0.1403, 0.0924, 0.2987, -0.003, 0.0823, -0.0094],
                    [-0.1768, 0.12, -0.2359, 0.0119, -0.024, 0.1089, -0.097, 0.0079, 0.0106, -0.3506],
                    [-0.1228, 0.0961, 0.139, 0.0072, -0.008, -0.0618, 0., 0.458, -0.1264, 0.0032],
                    [-0.1664, -0.2362, 0.0119, -0.0866, -0.211, 0.005, -0.0158, 0.051, 0.1973, 0.069],
                    [-0.1385, 0.0974, 0.0795, 0.3422, -0.1616, -0.1567, 0.1527, -0.0189, 0.0513, 0.0097],
                    [-0.1797, 0.1216, -0.2431, 0.0074, -0.0183, 0.1817, -0.2566, 0.0146, -0.0372, -0.0007],
                    [-0.165, -0.2268, 0.0072, 0.1019, 0.3011, -0.0281, -0.0335, 0.0244, 0.0911, 0.0322],
                    [-0.1647, -0.2338, 0.0113, -0.0493, -0.0814, 0.0024, 0.0206, -0.0692, -0.2795, -0.1135],
                    [-0.1513, 0.0767, -0.0498, 0.4283, -0.1305, 0.0737, 0.3194, -0.0151, 0.0811, 0.0685],
                    [-0.1292, -0.1716, 0.0159, -0.1305, -0.3642, 0.0657, 0.1177, -0.0863, -0.2945, 0.0722],
                    [-0.1546, 0.1158, 0.1941, -0.173, 0.0866, 0.2732, 0.1782, -0.1098, 0.0663, 0.0088],
                    [-0.1526, 0.1063, -0.2191, -0.2164, 0.0831, -0.2866, 0.2083, -0.0193, 0.0069, 0.256],
                    [-0.1546, 0.1225, 0.1995, -0.1765, 0.1059, 0.2562, 0.1741, 0.2687, -0.0187, -0.006],
                    [-0.1226, 0.0943, 0.1555, -0.1551, 0.0886, 0.2603, 0.2172, 0.3703, -0.0461, 0.0023],
                    [-0.1616, 0.1123, -0.2289, -0.1959, 0.0715, -0.2622, 0.2105, -0.0168, 0.019, 0.0623],
                    [-0.1615, 0.1285, 0.2084, -0.1584, 0.0923, 0.1926, 0.0896, 0.0544, 0.0182, -0.0086],
                    [-0.1771, 0.1217, -0.2443, -0.083, 0.0152, -0.0311, 0., -0.0016, 0.0157, -0.341],
                    [-0.1726, -0.2461, 0.0126, -0.1002, -0.233, 0.011, 0.0131, -0.0257, -0.1092, -0.057],
                    [-0.1743, 0.1381, 0.2159, -0.0326, 0.0077, -0.0662, -0.1282, -0.293, 0.05, -0.0053],
                    [-0.1739, 0.1365, 0.2003, 0.0859, -0.0619, -0.232, -0.1599, 0.2309, -0.1064, -0.0049],
                    [-0.1752, -0.2504, 0.0119, -0.0607, -0.0916, -0.0177, -0.0616, 0.0575, 0.1947, -0.0628],
                    [-0.1452, -0.2041, 0.0126, -0.1397, -0.3936, 0.0492, 0.0936, -0.0747, -0.2619, 0.0579],
                    [-0.1522, 0.1048, -0.2159, -0.1017, 0.0302, 0.0249, -0.1787, 0.0046, -0.049, 0.2499],
                    [-0.1775, 0.1403, 0.2121, 0.0554, -0.0479, -0.2257, -0.2011, 0.0832, -0.0765, -0.008],
                    [-0.1741, -0.2409, 0.0089, 0.0444, 0.1389, -0.0233, -0.0491, 0.0764, 0.2928, 0.0765]])

    np.set_printoptions(precision=6, linewidth=200)
    test(mat)

    r = main(mat)
    print(r[0])
    print(r[1])
    print(r[2])
