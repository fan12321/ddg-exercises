// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for (auto v: mesh.vertices()) {
        size_t index = v.getIndex();
        double dualArea = barycentricDualArea(v);
        tripletList.push_back(T(index, index, dualArea));
    }

    size_t nVertices = mesh.nVertices();
    SparseMatrix<double> resultMatrix(nVertices, nVertices);
    resultMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return resultMatrix; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for (auto e: mesh.edges()) {
        size_t index = e.getIndex();
        auto he1 = e.halfedge();
        auto he2 = he1.twin();
        double ratio = 0.5 * (cotan(he1) + cotan(he2));
        tripletList.push_back(T(index, index, ratio));
    }
    
    size_t nEdges = mesh.nEdges();
    SparseMatrix<double> resultMatrix(nEdges, nEdges);
    resultMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return resultMatrix; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for (auto f: mesh.faces()) {
        size_t index = f.getIndex();
        auto he = f.halfedge();
        Vertex a = he.tailVertex();
        Vertex b = he.tipVertex();
        Vertex c = he.next().tipVertex();

        Vector3 e1 = inputVertexPositions[b] - inputVertexPositions[a];
        Vector3 e2 = inputVertexPositions[c] - inputVertexPositions[a];

        double area = sqrt(
            (e1.y*e2.z - e1.z*e2.y) * (e1.y*e2.z - e1.z*e2.y) + 
            (e1.x*e2.z - e1.z*e2.x) * (e1.x*e2.z - e1.z*e2.x) + 
            (e1.x*e2.y - e1.y*e2.x) * (e1.x*e2.y - e1.y*e2.x)
        ) * 0.5;

        tripletList.push_back(T(index, index, 1.0 / area));
    }
    
    size_t nFaces = mesh.nFaces();
    SparseMatrix<double> resultMatrix(nFaces, nFaces);
    resultMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return resultMatrix; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for (auto e: mesh.edges()) {
        size_t edgeIndex = e.getIndex();
        size_t v0Index = e.firstVertex().getIndex();
        size_t v1Index = e.secondVertex().getIndex();
        tripletList.push_back(T(edgeIndex, v0Index, -1.0));
        tripletList.push_back(T(edgeIndex, v1Index, 1.0));
    }

    size_t nVertices = mesh.nVertices();
    size_t nEdges = mesh.nEdges();
    SparseMatrix<double> resultMatrix(nEdges, nVertices);
    resultMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return resultMatrix;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    for (auto f: mesh.faces()) {
        size_t faceIndex = f.getIndex();
        auto he = f.halfedge();
        size_t e0Index = he.edge().getIndex();
        size_t e1Index = he.next().edge().getIndex();
        size_t e2Index = he.next().next().edge().getIndex();
        double sign0 = (he.orientation())? 1.0 : -1.0;
        double sign1 = (he.next().orientation())? 1.0 : -1.0;
        double sign2 = (he.next().next().orientation())? 1.0 : -1.0;
        tripletList.push_back(T(faceIndex, e0Index, sign0));
        tripletList.push_back(T(faceIndex, e1Index, sign1));
        tripletList.push_back(T(faceIndex, e2Index, sign2));
    }

    size_t nEdges = mesh.nEdges();
    size_t nFaces = mesh.nFaces();
    SparseMatrix<double> resultMatrix(nFaces, nEdges);
    resultMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return resultMatrix;
}

} // namespace surface
} // namespace geometrycentral