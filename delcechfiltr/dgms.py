from itertools import combinations
import numpy as np
from scipy.spatial.distance import euclidean
from scipy.spatial import Delaunay
import delcechfiltr.tri
import gudhi

def simplices_and_cech_param(points):
    dim = len(points[0])
    if dim == 2:
        out_del = Delaunay(points)
        tri_del = out_del.simplices
        tri_del = sorted([sorted(tri) for tri in tri_del])
        edges_del = set()
        for v1, v2, v3 in tri_del:
            edges_del.add((v1, v2))
            edges_del.add((v1, v3))
            edges_del.add((v2, v3))
        edges_del = list(edges_del)

        param_edges = [euclidean(points[i], points[j]) / 2.0  for (i,j) in edges_del]
        param_tri = delcechfiltr.tri.cech_param_list(points, tri_del)

        sim0 = np.array([[v, -1, -1] for v in range(len(points))])
        sim1 = np.array([[ed[0], ed[1], -1] for ed in edges_del])
        simplices = np.vstack([sim0, sim1, tri_del])
        parameterization = np.hstack([np.zeros(len(points)),
                                     param_edges,
                                     param_tri])
        return simplices, parameterization
    elif dim == 3:
        out_del = Delaunay(points)
        tetra_del = out_del.simplices
        tetra_del = sorted([sorted(te) for te in tetra_del])
        edges_del = set()
        tri_del = set()
        for v1, v2, v3, v4 in tetra_del:
            edges_del.add((v1, v2))
            edges_del.add((v1, v3))
            edges_del.add((v1, v4))
            edges_del.add((v2, v3))
            edges_del.add((v2, v4))
            edges_del.add((v3, v4))

            tri_del.add((v1, v2, v3))
            tri_del.add((v1, v2, v4))
            tri_del.add((v1, v3, v4))
            tri_del.add((v2, v3, v4))

        edges_del = list(edges_del)
        tri_del = list(tri_del)

        param_edges = [euclidean(points[i], points[j]) / 2.0  for (i,j) in edges_del]
        param_tri = delcechfiltr.tri.cech_param_list(points, tri_del)
        param_tetra = delcechfiltr.tetra.cech_param_list(points, tetra_del)
        parameterization = np.hstack([np.zeros(len(points)),
                                     param_edges,
                                     param_tri,
                                     param_tetra])

        sim0 = np.array([[v, -1, -1, -1] for v in range(len(points))])
        sim1 = np.array([[ed[0], ed[1], -1, -1] for ed in edges_del])
        sim2 = np.array([[tri[0], tri[1], tri[2], -1]
                         for tri in tri_del])
        simplices = np.vstack([sim0, sim1, sim2, tetra_del])
        return simplices, parameterization
    else:
        print("elements of `points` must be 2 or 3 dimensional")
        return None

def cech(points, homology_coeff_field=2,
         min_persistence=0.000000001, persistence_dim_max=True):
    edges = list(combinations(range(len(points)), 2))
    tri = list(combinations(range(len(points)), 3))
    st = gudhi.SimplexTree()
    param_edges = [euclidean(points[i], points[j]) / 2.0  for (i,j) in edges]
    param_tri = delcechfiltr.tri.cech_param_list(points, tri)
    for i, v in enumerate(points):
        st.insert([i], 0.0)
    for e, r1 in zip(edges, param_edges):
        st.insert(e, r1)
    for t, r2 in zip(tri, param_tri):
        st.insert(t, r2)
    st.make_filtration_non_decreasing()
    dgms = st.persistence(homology_coeff_field=homology_coeff_field,
                          min_persistence=min_persistence,
                          persistence_dim_max=persistence_dim_max)
    dgm0 = np.array([p for dim, p in dgms if dim == 0 and p[1] != float("inf")])
    dgm0 = dgm0[np.argsort(dgm0[:,1])]
    dgm1 = np.array([p for dim, p in dgms if dim == 1])
    dgm1 = dgm1[np.argsort(dgm1[:,1])]
    return dgm0, dgm1

def delcech_2D(points, homology_coeff_field=2,
               min_persistence=0.000000001, persistence_dim_max=True):
    out_del = Delaunay(points)
    tri_del = out_del.simplices
    tri_del = sorted([sorted(tri) for tri in tri_del])
    edges_del = set()
    for v1, v2, v3 in tri_del:
        edges_del.add((v1, v2))
        edges_del.add((v1, v3))
        edges_del.add((v2, v3))
    edges_del = list(edges_del)

    st = gudhi.SimplexTree()
    param_edges = [euclidean(points[i], points[j]) / 2.0  for (i,j) in edges_del]
    param_tri = delcechfiltr.tri.cech_param_list(points, tri_del)
    for i, v in enumerate(points):
        st.insert([i], 0.0)
    for e, r1 in zip(edges_del, param_edges):
        st.insert(e, r1)
    for t, r2 in zip(tri_del, param_tri):
        st.insert(t, r2)
    st.make_filtration_non_decreasing()
    dgms = st.persistence(homology_coeff_field=homology_coeff_field,
                          min_persistence=min_persistence,
                          persistence_dim_max=persistence_dim_max)
    dgm0 = np.array([p for dim, p in dgms if dim == 0 and p[1] != float("inf")])
    dgm0 = dgm0[np.argsort(dgm0[:,1])]
    dgm1 = np.array([p for dim, p in dgms if dim == 1])
    dgm1 = dgm1[np.argsort(dgm1[:,1])]
    return dgm0, dgm1

def delcech_3D(points, homology_coeff_field=2,
               min_persistence=0.000000001, persistence_dim_max=True):
    out_del = Delaunay(points)
    tetra_del = out_del.simplices
    tetra_del = sorted([sorted(te) for te in tetra_del])
    edges_del = set()
    tri_del = set()
    for v1, v2, v3, v4 in tetra_del:
        edges_del.add((v1, v2))
        edges_del.add((v1, v3))
        edges_del.add((v1, v4))
        edges_del.add((v2, v3))
        edges_del.add((v2, v4))
        edges_del.add((v3, v4))

        tri_del.add((v1, v2, v3))
        tri_del.add((v1, v2, v4))
        tri_del.add((v1, v3, v4))
        tri_del.add((v2, v3, v4))

    edges_del = list(edges_del)
    tri_del = list(tri_del)

    st = gudhi.SimplexTree()
    param_edges = [euclidean(points[i], points[j]) / 2.0  for (i,j) in edges_del]
    param_tri = delcechfiltr.tri.cech_param_list(points, tri_del)
    for i, v in enumerate(points):
        st.insert([i], 0.0)
    for e, r in zip(edges_del, param_edges):
        st.insert(e, r)
    for t, r in zip(tri_del, param_tri):
        st.insert(t, r)
    st.make_filtration_non_decreasing()
    dgms = st.persistence(homology_coeff_field=homology_coeff_field,
                          min_persistence=min_persistence,
                          persistence_dim_max=persistence_dim_max)
    dgm0 = np.array([p for dim, p in dgms if dim == 0 and p[1] != float("inf")])
    dgm0 = dgm0[np.argsort(dgm0[:,1])]
    dgm1 = np.array([p for dim, p in dgms if dim == 1])
    dgm1 = dgm1[np.argsort(dgm1[:,1])]
    return dgm0, dgm1
