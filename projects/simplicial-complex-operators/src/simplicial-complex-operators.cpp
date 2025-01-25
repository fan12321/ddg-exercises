// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    typedef Eigen::Triplet<float> T;
    std::vector<T> tripletList;
    size_t nRows = mesh->nEdges();
    size_t nCols = mesh->nVertices();

    tripletList.reserve(nRows * 2);
    for (Edge e : mesh->edges()) {
        auto edgeIdx = geometry->edgeIndices[e];
        auto vertexIdx_1 = geometry->vertexIndices[e.firstVertex()];
        auto vertexIdx_2 = geometry->vertexIndices[e.secondVertex()];
        tripletList.push_back(T(edgeIdx, vertexIdx_1, 1));
        tripletList.push_back(T(edgeIdx, vertexIdx_2, 1));
    }
    SparseMatrix<size_t> result(nRows, nCols);
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    typedef Eigen::Triplet<float> T;
    std::vector<T> tripletList;
    size_t nRows = mesh->nFaces();
    size_t nCols = mesh->nEdges();

    tripletList.reserve(nRows * 3);
    for (Face f : mesh->faces()) {
        auto faceIdx = geometry->faceIndices[f];
        Halfedge starting_HE = f.halfedge();
        Halfedge current_HE = starting_HE;
        do {
            auto edgeIdx = geometry->edgeIndices[current_HE.edge()];
            tripletList.push_back(T(faceIdx, edgeIdx, 1));
            current_HE = current_HE.next();
        } while (current_HE != starting_HE);
    }
    SparseMatrix<size_t> result(nRows, nCols);
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    size_t length = A0.cols();
    Vector<size_t> result(length);
    for (auto i=0; i<length; i++) {
        result[i] = subset.vertices.count(i);
    }
    return result;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    size_t length = A1.cols();
    Vector<size_t> result(length);
    for (auto i=0; i<length; i++) {
        result[i] = subset.edges.count(i);
    }
    return result;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    size_t length = A1.rows();
    Vector<size_t> result(length);
    for (auto i=0; i<length; i++) {
        result[i] = subset.faces.count(i);
    }
    return result;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    MeshSubset result = subset.deepCopy();

    for (auto edgeIdx: subset.edges) {
        size_t startingIdx_A1 = A1.outerIndexPtr()[edgeIdx];
        size_t endingIdx_A1 = A1.outerIndexPtr()[edgeIdx+1];
        for (size_t j=startingIdx_A1; j<endingIdx_A1; j++) {
            size_t faceIdx = A1.innerIndexPtr()[j];
            result.addFace(faceIdx);
        }
    }

    for (auto vertexIdx: subset.vertices) {
        // add edges containing the vertex to the result
        size_t startingIdx_A0 = A0.outerIndexPtr()[vertexIdx];
        size_t endingIdx_A0 = A0.outerIndexPtr()[vertexIdx+1];
        for (size_t i=startingIdx_A0; i<endingIdx_A0; i++) {
            size_t edgeIdx = A0.innerIndexPtr()[i];
            // skip if the edge is in the original subset since it has already been checked
            if (subset.edges.count(edgeIdx)) continue;
            
            result.addEdge(edgeIdx);

            // add faces containing the edge (hence containing the vertex) to the result
            size_t startingIdx_A1 = A1.outerIndexPtr()[edgeIdx];
            size_t endingIdx_A1 = A1.outerIndexPtr()[edgeIdx+1];
            for (size_t j=startingIdx_A1; j<endingIdx_A1; j++) {
                auto faceIdx = A1.innerIndexPtr()[j];
                result.addFace(faceIdx);
            }
        }
    }
    return result; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset result = subset.deepCopy();

    // face
    std::vector<size_t> A1ColumnIndex(A1.nonZeros());
    std::vector<size_t> A0ColumnIndex(A0.nonZeros());
    // A1 column indices
    size_t outerIndex = 0;
    for (size_t innerIndex=0; innerIndex<A1.nonZeros(); innerIndex++) {
        if (innerIndex == A1.outerIndexPtr()[outerIndex+1]) outerIndex += 1;
        A1ColumnIndex[innerIndex] = outerIndex;
    }
    // A0 column indices
    outerIndex = 0;
    for (size_t innerIndex=0; innerIndex<A0.nonZeros(); innerIndex++) {
        if (innerIndex == A0.outerIndexPtr()[outerIndex+1]) outerIndex += 1;
        A0ColumnIndex[innerIndex] = outerIndex;
    }

    for (size_t faceEdgeIndex=0; faceEdgeIndex<A1.nonZeros(); faceEdgeIndex++) {
        size_t faceIdx = A1.innerIndexPtr()[faceEdgeIndex];
        size_t edgeIdx = A1ColumnIndex[faceEdgeIndex];
        if (subset.faces.count(faceIdx)) {
            result.addEdge(edgeIdx);
            for (size_t edgeVertexIndex=0; edgeVertexIndex<A0.nonZeros(); edgeVertexIndex++) {
                if (A0.innerIndexPtr()[edgeVertexIndex] == edgeIdx) result.addVertex(A0ColumnIndex[edgeVertexIndex]);
            }
        }
    }

    // edge
    for (size_t edgeIdx: subset.edges) {
        for (size_t edgeVertexIndex=0; edgeVertexIndex<A0.nonZeros(); edgeVertexIndex++) {
            if (A0.innerIndexPtr()[edgeVertexIndex] == edgeIdx) result.addVertex(A0ColumnIndex[edgeVertexIndex]);
        }
    }

    return result; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    // closure of star minus star of closure
    auto starOfClosure = star(closure(subset));
    auto result = closure(star(subset));
    result.deleteSubset(starOfClosure);
    return result; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    // face
    std::vector<size_t> A1ColumnIndex(A1.nonZeros());
    std::vector<size_t> A0ColumnIndex(A0.nonZeros());
    // A1 column indices
    size_t outerIndex = 0;
    for (size_t innerIndex=0; innerIndex<A1.nonZeros(); innerIndex++) {
        if (innerIndex == A1.outerIndexPtr()[outerIndex+1]) outerIndex += 1;
        A1ColumnIndex[innerIndex] = outerIndex;
    }
    // A0 column indices
    outerIndex = 0;
    for (size_t innerIndex=0; innerIndex<A0.nonZeros(); innerIndex++) {
        if (innerIndex == A0.outerIndexPtr()[outerIndex+1]) outerIndex += 1;
        A0ColumnIndex[innerIndex] = outerIndex;
    }

    // for each face, the edges should also be in the subset
    for (size_t faceEdgeIndex=0; faceEdgeIndex<A1.nonZeros(); faceEdgeIndex++) {
        size_t faceIdx = A1.innerIndexPtr()[faceEdgeIndex];
        size_t edgeIdx = A1ColumnIndex[faceEdgeIndex];
        if (subset.faces.count(faceIdx) && !subset.edges.count(edgeIdx)) return false;
    }

    // for each edge, the vertices should also be in the subset
    for (size_t edgeVertexIndex=0; edgeVertexIndex<A0.nonZeros(); edgeVertexIndex++) {
        size_t edgeIdx = A0.innerIndexPtr()[edgeVertexIndex];
        size_t vertexIdx = A0ColumnIndex[edgeVertexIndex];
        if (subset.edges.count(edgeIdx) && !subset.vertices.count(vertexIdx)) return false; 
    }

    return true; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (!isComplex(subset)) return -1;
    else if (
        subset.vertices.empty() && 
        subset.edges.empty() && 
        subset.faces.empty()
    ) return -1;
    else if (
        !subset.vertices.empty() &&
        subset.edges.empty() && 
        subset.faces.empty()
    ) return 0;
    else if (
        !subset.vertices.empty() &&
        !subset.edges.empty() && 
        subset.faces.empty()
    ) {
        // check if each 0-simplex contained in an edge
        for (auto vertexIdx: subset.vertices) {
            auto startingIdx_A0 = A0.outerIndexPtr()[vertexIdx];
            auto endingIdx_A0 = A0.outerIndexPtr()[vertexIdx+1];
            bool containedByEdge = false;
            for (size_t i=startingIdx_A0; i<endingIdx_A0; i++) {
                auto edgeIdx = A0.innerIndexPtr()[i];
                if (subset.edges.count(edgeIdx)) {
                    containedByEdge = true;
                    break;
                }
            }
            if (!containedByEdge) return -1;
        }
        return 1;
    }
    else {
        // check if each 0-simplex contained in an edge
        for (auto vertexIdx: subset.vertices) {
            auto startingIdx_A0 = A0.outerIndexPtr()[vertexIdx];
            auto endingIdx_A0 = A0.outerIndexPtr()[vertexIdx+1];
            bool containedByEdge = false;
            for (size_t i=startingIdx_A0; i<endingIdx_A0; i++) {
                auto edgeIdx = A0.innerIndexPtr()[i];
                if (subset.edges.count(edgeIdx)) {
                    containedByEdge = true;
                    break;
                }
            }
            if (!containedByEdge) return -1;
        }

        // check if each 1-simplex contained
        for (auto edgeIdx: subset.edges) {
            auto startingIdx_A1 = A1.outerIndexPtr()[edgeIdx];
            auto endingIdx_A1 = A1.outerIndexPtr()[edgeIdx+1];
            bool containedByFace = false;
            for (size_t i=startingIdx_A1; i<endingIdx_A1; i++) {
                auto faceIdx = A1.innerIndexPtr()[i];
                if (subset.faces.count(faceIdx)) {
                    containedByFace = true;
                    break;
                }
            }
            if (!containedByFace) return -1;
        }
        return 2;
    }
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    int pureComplexDegree = isPureComplex(subset);
    switch (pureComplexDegree) {
    case -1: {
        // not defined
        return MeshSubset();
    }
    case 0: {
        // null set
        return MeshSubset();
    }
    case 1: {
        // vertices that is only contained by one edge
        MeshSubset result({}, {}, {});
        for (auto vertexIdx: subset.vertices) {
            auto startingIdx_A0 = A0.outerIndexPtr()[vertexIdx];
            auto endingIdx_A0 = A0.outerIndexPtr()[vertexIdx+1];
            int numberOfConnectingEdges = 0;
            for (auto i=startingIdx_A0; i<endingIdx_A0; i++) {
                auto edgeIdx = A0.innerIndexPtr()[i];
                if (subset.edges.count(edgeIdx)) {
                    numberOfConnectingEdges += 1;
                }
            }
            if (numberOfConnectingEdges == 1) result.addVertex(vertexIdx);
        }
        return result;
    }
    case 2: {
        MeshSubset result({}, {}, {});
        for (auto edgeIdx: subset.edges) {
            auto startingIdx_A1 = A1.outerIndexPtr()[edgeIdx];
            auto endingIdx_A1 = A1.outerIndexPtr()[edgeIdx+1];
            int numberOfConnectingFaces = 0;
            for (auto i=startingIdx_A1; i<endingIdx_A1; i++) {
                auto faceIdx = A1.innerIndexPtr()[i];
                if (subset.faces.count(faceIdx)) {
                    numberOfConnectingFaces += 1;
                }
            }
            if (numberOfConnectingFaces == 1) result.addEdge(edgeIdx);
        }
        result = closure(result);
        return result;
    }
    default: {
        return MeshSubset();
    }
    }
}