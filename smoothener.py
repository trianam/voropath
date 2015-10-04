import numpy as np

class Smoothener:
    def __init__(self):
        pass

    def plot(self, V, plotter):
        k = 3
        den = float(len(V) - k + 1)
        p = []
        for i in range(k):
            p[:0] = [1.]
        for i in range(len(V)-k)[::-1]:
            p[:0] = [(float(i)+1.) / den]
        for i in range(k):
            p[:0] = [0.]

        p = np.array(p)
        plotter.plot(V[:,0], V[:,1], V[:,2], 'r--')
        spline = self._bSpline(V, k, p)
        plotter.plot(spline[:,0], spline[:,1], spline[:,2], 'r', lw=2)

    def _bSplineBaseFactory(self,part,k):
        """
        Define the B-splines base for the given extended partition.

        Arguments:
        - part: the extended partition of the parametric domain.
        - k: the order

        Returns: the function that give a base for a given i and h.
        """
        if (part[len(part)-k+1 : len(part)] == part[len(part)-k : len(part)-1]).all():
            iEnd = len(part)-k
        else:
            iEnd = len(part) -1

        def N(i,h=k):
            """
            Define a single base of B-splines for a given i and h.

            Arguments:
            - i: the index of the base
            - h: the order of the base

            Returns:
            - the function in t of the base
            """
            def bSpline(t):
                if (h > 1):
                    n1 = N(i,h-1)
                    n2 = N(i+1, h-1)

                    a1 = ((t - part[i]) * n1(t))
                    a2 = (part[i+h-1] - part[i])
                    if(a1 == 0.):
                        a = 0.
                    else:
                        a = (a1 / a2)

                    b1 = ((part[i+h] - t) * n2(t))
                    b2 = (part[i+h] - part[i+1])
                    if(b1 == 0.):
                        b = 0.
                    else:
                        b = (b1 / b2)

                    return a + b
                else:
                    #if ((part[i] <= t) and (t < part[i+1])):
                    #if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == len(part)-k) and (t == part[i+1])))):
                    #if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == len(part)-1) and (t == part[i+1])))):
                    if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == iEnd) and (t == part[iEnd])))):

                        return 1.
                    else:
                        return 0.
            return bSpline
        return N

    def _bSplineFun(self, V, k, dom=(0.,1.), part=None, m=None):
        """
        Define a B-spline function using the given control vectors, of
        the given order, and using either the given extended partition
        or creating one uniform and clamped in the given domain and
        molteplicity.

        Arguments:
        - V: the control points vector;
        - k: the order;
        - dom: the parametric domain, mandatory if a partition is not
        given, but with default value (0,1);
        - part: the extended partition of the parametric domain,
        mandatory if the multiplicity and domain are not given;
        - m: the multiplicity, mandatory if the partition is not given.

        Returns: the function in t of the curve
        """
        if (part == None):
            part = np.concatenate((np.zeros(k)+dom[0], np.repeat(np.linspace(dom[0], dom[1], len(m)+2)[1:-1], m), np.zeros(k)+dom[1]))
        else:
            dom = (part[k-1], part[-k])

        N = self._bSplineBaseFactory(part,k)
        return lambda t: sum(V[i]*N(i)(t) for i in range(len(V)))

    def _bSpline(self, V, k, part, step = 0.01):
        """
        calculate the points of the B-spline curve defined by the
        given control points, the given extended partition, and of
        the given order.

        Arguments:
        - V: the control points vector;
        - k: the order;
        - part: the extended partition of the parametric domain;
        - step: the steps in the domain for creating the points.

        Returns:
        the array of points in the space of the curve.
        """
        C = self._bSplineFun(V,k,part=part)

        P = np.empty((0,V.shape[1]), dtype=float);
        for t in np.arange(part[k-1],part[-k]+step,step):
            P = np.append(P, [C(t)], axis=0)

        return P
