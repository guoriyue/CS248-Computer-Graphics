#include <queue>
#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include "../geometry/halfedge.h"
#include "debug.h"

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementaiton, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {
    // (void)v;
    // return std::nullopt;
    if(v->on_boundary()) {
        (void)v;
        return std::nullopt;
    }

    std::vector<Halfedge_Mesh::HalfedgeRef> halfedges_to_erase;
    Halfedge_Mesh::HalfedgeRef h = v->halfedge();
    halfedges_to_erase.push_back(h);
    
    while (h->twin()->next() != v->halfedge())
    {
        h = h->twin()->next();
        halfedges_to_erase.push_back(h);
    }

    h = v->halfedge()->next();
    
    erase(v);
    auto f = Halfedge_Mesh::erase_edge(halfedges_to_erase[0]->edge());
    for(size_t i = 1; i < halfedges_to_erase.size(); i++) {
        f = Halfedge_Mesh::erase_edge(halfedges_to_erase[i]->edge());
    }
    
    return f;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {
    // (void)e;
    // return std::nullopt;
    if(e->on_boundary()) {
        (void)e;
        return std::nullopt;
    }
    
    
    Halfedge_Mesh::FaceRef face0 = e->halfedge()->face();
    Halfedge_Mesh::FaceRef face1 = e->halfedge()->twin()->face();
    Halfedge_Mesh::HalfedgeRef halfedge0 = e->halfedge();
    Halfedge_Mesh::HalfedgeRef halfedge1 = e->halfedge()->twin();

    
    if (e->halfedge()->face() == e->halfedge()->twin()->face())
    {
        // only one edge
        // both pair of next and next_prev is valid, halfedge and halfedge twin is a loop
        if (e->halfedge()->next() == e->halfedge()->twin() || e->halfedge()->twin()->next() == e->halfedge())
        {
            Halfedge_Mesh::HalfedgeRef halfedge0_next = e->halfedge()->next();
            Halfedge_Mesh::HalfedgeRef halfedge0_next_prev_prev = e->halfedge()->next()->next();
            while (halfedge0_next_prev_prev->next()->next() != halfedge0_next) {
                halfedge0_next_prev_prev = halfedge0_next_prev_prev->next();
            }
            Halfedge_Mesh::HalfedgeRef halfedge1_next = e->halfedge()->twin()->next();
            Halfedge_Mesh::HalfedgeRef halfedge1_next_prev_prev = e->halfedge()->twin()->next()->next();
            while (halfedge1_next_prev_prev->next()->next() != halfedge1_next) {
                halfedge1_next_prev_prev = halfedge1_next_prev_prev->next();
            }
            if (e->halfedge()->next() == e->halfedge()->twin())
            {
                halfedge0_next_prev_prev->next() = halfedge0_next->next();
                e->halfedge()->face()->halfedge() = halfedge0_next_prev_prev;
                halfedge0_next_prev_prev->next()->vertex()->halfedge() = halfedge0_next->next();
            }
            else
            {
                halfedge1_next_prev_prev->next() = halfedge1_next->next();
                e->halfedge()->face()->halfedge() = halfedge1_next_prev_prev;
                halfedge1_next_prev_prev->next()->vertex()->halfedge() = halfedge1_next->next();
            }
            
            // erase edge, halfedge, face
            erase(e);
            erase(halfedge0);
            erase(halfedge1);
            return face0;
        }
    }

    // erase halfedge0 and halfedge1, make the vertex of halfedge0 point to halfedge1->next()
    halfedge0->vertex()->halfedge() = halfedge1->next();
    halfedge1->vertex()->halfedge() = halfedge0->next();

    face0->halfedge() = halfedge0->next();

    // h0->prev->next = h1->next
    Halfedge_Mesh::HalfedgeRef h = halfedge0->next();
    // h->face() = face0;
    while (h->next() != halfedge0) {
        h = h->next();
    } 
    h->next() = halfedge1->next();
    // h1->prev->next = h0->next
    h = halfedge1->next();
    while (h->next() != halfedge1) {
        h->face() = face0;
        h = h->next();
    }
    h->face() = face0;
    h->next() = halfedge0->next();

    // erase edge, halfedge, face
    erase(e);
    erase(halfedge0);
    erase(halfedge1);
    // erase face1, when there are only 1 face left, do not erase it
    if (face1 != face0)
        erase(halfedge1->face());
    return face0;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {
    bool boundary= e->on_boundary();
    
    if (boundary) {
        (void)e;
        return std::nullopt;
    }
    
    if (e->halfedge()->face() == e->halfedge()->twin()->face())
    {
        (void)e;
        return std::nullopt;
    }
    // multiple edges error
    // we only handle two triangles here, so if there are more than two triangles, return nullopt
    // more than two triangles: a different face from the ones we are probably going to erase
    std::unordered_set<VertexRef> vertices;
    Halfedge_Mesh::HalfedgeRef halfedge0 = e->halfedge();
    Halfedge_Mesh::HalfedgeRef halfedge1 = e->halfedge()->twin();
    Halfedge_Mesh::FaceRef face0 = e->halfedge()->face();
    Halfedge_Mesh::FaceRef face1 = e->halfedge()->twin()->face();
    Halfedge_Mesh::HalfedgeRef h = halfedge0;
    
    while (h->twin()->next() != halfedge0) {
        // add vertices to the set
        // ignore vertices of the edge we are collapsing
        if (h != halfedge0)
            vertices.insert(h->twin()->vertex());
        vertices.insert(h->twin()->vertex());
        h = h->twin()->next();
    }

    h = halfedge1;
    while (h->twin()->next() != halfedge1) {
        if (vertices.find(h->twin()->vertex()) != vertices.end() && (h->face() != face0 && h->face() != face1) && (h->twin()->face() != face0 && h->twin()->face() != face1))
        {
            // check if there are duplicate vertices in the set
            // need to use twin vertex to get the vertex of the other end of the edge
            (void)e;
            return std::nullopt;
        }
        if (h != halfedge1)
            vertices.insert(h->twin()->vertex());
        h = h->twin()->next();
    }
    
    // Halfedge_Mesh::HalfedgeRef halfedge0 = e->halfedge();
    // Halfedge_Mesh::HalfedgeRef halfedge1 = e->halfedge()->twin();
    halfedge0 = e->halfedge();
    halfedge1 = e->halfedge()->twin();
    if (e->halfedge()->is_boundary())
    {
        halfedge0 = e->halfedge()->twin();
    }
    // if either of the polygons containing e was a triangle, it will be replaced by an edge
    Halfedge_Mesh::VertexRef vertex0 = halfedge0->vertex();
    Halfedge_Mesh::VertexRef vertex1 = halfedge1->vertex();
    Halfedge_Mesh::VertexRef mid = new_vertex();
    mid->pos = e->center();


    
    Halfedge_Mesh::HalfedgeRef halfedge0_next = halfedge0->next();
    Halfedge_Mesh::HalfedgeRef halfedge1_next = halfedge1->next();
    Halfedge_Mesh::HalfedgeRef halfedge0_prev = halfedge0->next();
    Halfedge_Mesh::HalfedgeRef halfedge1_prev = halfedge1->next();
    while (halfedge0_prev->next() != halfedge0)
    {
        halfedge0_prev = halfedge0_prev->next();
    }
    while (halfedge1_prev->next() != halfedge1)
    {
        halfedge1_prev = halfedge1_prev->next();
    }

    if (!boundary)
    {
        // update the halfedge of the face because these two halfedges will collapse
        halfedge0->face()->halfedge() = halfedge0_next;
        halfedge1->face()->halfedge() = halfedge1_next;

        // set new vertex, perform the collapse
        Halfedge_Mesh::HalfedgeRef h = halfedge0_next;
        while (h->twin()->next() != halfedge0_next)
        {
            h->vertex() = mid;
            h = h->twin()->next();
        }
        h->vertex() = mid;

        h = halfedge1_next;
        while (h->twin()->next() != halfedge1_next)
        {
            h->vertex() = mid;
            h = h->twin()->next();
        }
        h->vertex() = mid;

        halfedge0_prev->next() = halfedge0_next;
        halfedge1_prev->next() = halfedge1_next;

        // erase the old vertex, edge, halfedge, face
        erase(vertex0);
        erase(vertex1);
        erase(e);
        erase(halfedge0);
        erase(halfedge1);

        // for triangle, collapse to an edge
        if(halfedge0_next->next() == halfedge0_prev)
        {
            // erase the halfedges inside the triangle, keep prev and erase next
            halfedge0_prev->twin()->twin() = halfedge0_next->twin();
            halfedge0_next->twin()->twin() = halfedge0_prev->twin();
            
            halfedge0_prev->twin()->edge()->halfedge() = halfedge0_prev->twin();
            halfedge0_prev->edge()->halfedge() = halfedge0_next->twin();
            halfedge0_prev->twin()->vertex()->halfedge() = halfedge0_prev->twin();
            halfedge0_prev->vertex()->halfedge() = halfedge0_next->twin();

            halfedge0_next->twin()->edge() = halfedge0_prev->twin()->edge();
            halfedge0_next->twin()->edge()->halfedge() = halfedge0_next->twin();

            mid->halfedge() = halfedge0_prev->twin();
            erase(halfedge0_next->face());
            erase(halfedge0_next->edge());
            erase(halfedge0_next);
            erase(halfedge0_prev);
            // a live halfedge's twin was erased, so we need to update the twin
        }
        if(halfedge1_next->next() == halfedge1_prev)
        {
            halfedge1_prev->twin()->twin() = halfedge1_next->twin();
            halfedge1_next->twin()->twin() = halfedge1_prev->twin();
            
            halfedge1_prev->twin()->edge()->halfedge() = halfedge1_prev->twin();
            halfedge1_prev->edge()->halfedge() = halfedge1_next->twin();
            halfedge1_prev->twin()->vertex()->halfedge() = halfedge1_prev->twin();
            halfedge1_prev->vertex()->halfedge() = halfedge1_next->twin();

            halfedge1_next->twin()->edge() = halfedge1_prev->twin()->edge();
            halfedge1_next->twin()->edge()->halfedge() = halfedge1_next->twin();

            mid->halfedge() = halfedge1_prev->twin();
            erase(halfedge1_next->face());
            erase(halfedge1_next->edge());
            erase(halfedge1_next);
            erase(halfedge1_prev);
        }
        return mid;
    }
    else
    {
        // advanced task: collapse a boundary edge

    }
    return std::nullopt;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {

    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {

    // when the edge is on the boundary, do nothing.
    if(e->on_boundary()) return std::nullopt;

    /*
    The following code only works for triangle mesh
    // PHASE I: Collect elements
    // Collect all halfedges
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h1 = h0->next();
    HalfedgeRef h2 = h1->next();
    HalfedgeRef h3 = h0->twin();
    HalfedgeRef h4 = h3->next();
    HalfedgeRef h5 = h4->next();
    HalfedgeRef h6 = h1->twin();
    HalfedgeRef h7 = h2->twin();
    HalfedgeRef h8 = h4->twin();
    HalfedgeRef h9 = h5->twin();

    // Collect all vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    VertexRef v2 = h5->vertex();
    VertexRef v3 = h2->vertex();

    // Collect all edges
    EdgeRef e1 = h5->edge();
    EdgeRef e2 = h4->edge();
    EdgeRef e3 = h2->edge();
    EdgeRef e4 = h1->edge();

    // Collect all faces
    FaceRef f0 = h0->face();
    FaceRef f1 = h3->face();

    // PHASE II: Allocate new elements
    // Nothing for flip edge

    // PHASE III: Reassign Elements
    // Reassign halfedges

    h0->next() = h1;
    h0->twin() = h3;
    h0->vertex() = v2;
    h0->edge() = e;
    h0->face() = f0;

    h1->next() = h2;
    h1->twin() = h7;
    h1->vertex() = v3;
    h1->edge() = e3;
    h1->face() = f0;

    h2->next() = h0;
    h2->twin() = h8;
    h2->vertex() = v0;
    h2->edge() = e2;
    h2->face() = f0;

    h3->next() = h4;
    h3->twin() = h0;
    h3->vertex() = v3;
    h3->edge() = e;
    h3->face() = f1;

    h4->next() = h5;
    h4->twin() = h9;
    h4->vertex() = v2;
    h4->edge() = e1;
    h4->face() = f1;

    h5->next() = h3;
    h5->twin() = h6;
    h5->vertex() = v1;
    h5->edge() = e4;
    h5->face() = f1;

    h6->next() = h6->next();
    h6->twin() = h5;
    h6->vertex() = v3;
    h6->edge() = e4;
    h6->face() = h6->face();

    h7->next() = h7->next();
    h7->twin() = h1;
    h7->vertex() = v0;
    h7->edge() = e3;
    h7->face() = h7->face();

    h8->next() = h8->next();
    h8->twin() = h2;
    h8->vertex() = v2;
    h8->edge() = e2;
    h8->face() = h8->face();

    h9->next() = h9->next();
    h9->twin() = h4;
    h9->vertex() = v1;
    h9->edge() = e1;
    h9->face() = h9->face();

    // Reassign vertices
    v0->halfedge() = h2;
    v1->halfedge() = h5;
    v2->halfedge() = h4;
    v3->halfedge() = h3;

    // Reassign edges
    e->halfedge() = h0;
    e1->halfedge() = h4;
    e2->halfedge() = h2;
    e3->halfedge() = h1;
    e4->halfedge() = h5;

    // Reassign faces
    f0->halfedge() = h0;
    f1->halfedge() = h3;

    // PHASE IV: Delete unused elements
    // Do nothing
    */

    // PHASE I: Collect elements
    // Collect halfedges
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h1 = h0->next();
    HalfedgeRef h2 = h1->next();
    HalfedgeRef h3 = h0;
    while(true) {
        h3 = h3->next();
        if(h3->next() == h0) break;
    }

    HalfedgeRef h4 = h0->twin();
    HalfedgeRef h5 = h4->next();
    HalfedgeRef h6 = h5->next();
    HalfedgeRef h7 = h4;
    while(true) {
        h7 = h7->next();
        if(h7->next() == h4) break;
    }
    // Collect vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h1->vertex();
    VertexRef v2 = h2->vertex();
    VertexRef v3 = h6->vertex();

    // Collect faces
    FaceRef f0 = h0->face();
    FaceRef f1 = h4->face();

    // PHASE III: Reassign Elements
    h0->next() = h2;
    h0->vertex() = v3;
    h1->next() = h4;
    h1->face() = f1;
    h3->next() = h5;

    h4->next() = h6;
    h4->vertex() = v2;
    h5->next() = h0;
    h5->face() = f0;
    h7->next() = h1;

    v0->halfedge() = h5;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h6;
    
    e->halfedge() = h0;

    f0->halfedge() = h0;
    f1->halfedge() = h4;

    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    // (void)e;
    // return std::nullopt;

    // this method is for triangle meshes only!
    if(e->halfedge()->face()->degree() != 3 || e->halfedge()->twin()->face()->degree() != 3) {
        (void)e;
        return std::nullopt;
    }

    // allocate space for new edges, halfedges, and faces
    std::vector<Halfedge_Mesh::EdgeRef> edges;
    std::vector<Halfedge_Mesh::HalfedgeRef> halfedges;
    std::vector<Halfedge_Mesh::FaceRef> faces;
    Halfedge_Mesh::VertexRef mid = new_vertex();
    mid->is_new = true;
    mid->pos = e->center();

    // halfedges we need to update
    Halfedge_Mesh::HalfedgeRef h0 = e->halfedge();
    Halfedge_Mesh::HalfedgeRef h1 = h0->next();
    Halfedge_Mesh::HalfedgeRef h2 = h1->next();
    Halfedge_Mesh::HalfedgeRef h3 = h0->twin();
    Halfedge_Mesh::HalfedgeRef h4 = h3->next();
    Halfedge_Mesh::HalfedgeRef h5 = h4->next();

    // halfedges for mid vertex and two old vertices
    mid->halfedge() = h0;
    h0->vertex()->halfedge() = h4;
    h3->vertex()->halfedge() = h1;
    h1->vertex()->halfedge() = h1;
    h2->vertex()->halfedge() = h2;
    h4->vertex()->halfedge() = h4;
    h5->vertex()->halfedge() = h5;

    // halfedges for faces
    faces.push_back(new_face());
    faces.push_back(new_face());
    
    // for faces
    faces[0]->halfedge() = h4;
    faces[1]->halfedge() = h2;
    h0->face()->halfedge() = h0;
    h3->face()->halfedge() = h3;

    // faces for other edges
    h1->face() = h0->face();
    h2->face() = faces[1];
    h4->face() = faces[0];
    h5->face() = h3->face();

    // assign new edges and halfedges
    EdgeRef e_new = new_edge();
    e_new->is_new = true;
    edges.push_back(new_edge());
    edges.push_back(new_edge());
    edges.push_back(e_new);
    halfedges.push_back(new_halfedge());
    halfedges.push_back(new_halfedge());
    halfedges.push_back(new_halfedge());
    halfedges.push_back(new_halfedge());
    halfedges.push_back(new_halfedge());
    halfedges.push_back(new_halfedge());
    edges[0]->halfedge() = halfedges[0];
    edges[0]->halfedge()->twin() = halfedges[1];
    edges[1]->halfedge() = halfedges[2];
    edges[1]->halfedge()->twin() = halfedges[3];
    edges[2]->halfedge() = halfedges[4];
    edges[2]->halfedge()->twin() = halfedges[5];


    // 8 halfedges
    h0->Halfedge_Mesh::Halfedge::set_neighbors(h1, h3, mid, e, h0->face());
    h3->Halfedge_Mesh::Halfedge::set_neighbors(h4, h0, h3->vertex(), e, h3->face());
    halfedges[0]->Halfedge_Mesh::Halfedge::set_neighbors(h5, halfedges[1], mid, edges[0], h3->face());
    halfedges[1]->Halfedge_Mesh::Halfedge::set_neighbors(halfedges[2], halfedges[0], h5->vertex(), edges[0], faces[0]);
    halfedges[2]->Halfedge_Mesh::Halfedge::set_neighbors(h4, halfedges[3], mid, edges[1], faces[0]);
    halfedges[3]->Halfedge_Mesh::Halfedge::set_neighbors(halfedges[4], halfedges[2], h4->vertex(), edges[1], faces[1]);
    halfedges[4]->Halfedge_Mesh::Halfedge::set_neighbors(h2, halfedges[5], mid, edges[2], faces[1]);
    halfedges[5]->Halfedge_Mesh::Halfedge::set_neighbors(h0, halfedges[4], h2->vertex(), edges[2], h0->face());
    h0->next() = h1;
    h1->next() = halfedges[5];
    h2->next() = halfedges[3];
    h3->next() = halfedges[0];
    h4->next() = halfedges[1];
    h5->next() = h3;

    return mid;
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)v;
    return std::nullopt;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    std::vector<HalfedgeRef> orignal_halfedges;
    std::vector<VertexRef> orignal_vertices;

    std::vector<EdgeRef> new_inward_edges;
    std::vector<EdgeRef> new_face_edges;

    std::vector<HalfedgeRef> inward_down_halfedges;
    std::vector<HalfedgeRef> inward_up_halfedges;

    std::vector<HalfedgeRef> face_cw_halfedges;
    std::vector<HalfedgeRef> face_ccw_halfedges;

    std::vector<VertexRef> new_vertices;
    std::vector<FaceRef> new_faces;

    HalfedgeRef he = f->halfedge();
    HalfedgeRef he0 = he;
    while(true) {
        orignal_halfedges.push_back(he);
        VertexRef v = he->vertex();
        orignal_vertices.push_back(v);
        he = he->next();
        if(he == he0) break;
    }

    int degree = f->degree();
    for(int i = 0; i < degree; i++) {

        inward_down_halfedges.push_back(new_halfedge());
        inward_up_halfedges.push_back(new_halfedge());
        face_cw_halfedges.push_back(new_halfedge());
        face_ccw_halfedges.push_back(new_halfedge());

        new_inward_edges.push_back(new_edge());
        new_face_edges.push_back(new_edge());
        new_vertices.push_back(new_vertex());
        new_faces.push_back(new_face());
    }

    for(int i = 0; i < degree; i++) {
        inward_down_halfedges[i]->next() = face_cw_halfedges[i == 0 ? degree - 1 : i - 1];
        inward_down_halfedges[i]->twin() = inward_up_halfedges[i];
        inward_down_halfedges[i]->vertex() = orignal_vertices[i];
        inward_down_halfedges[i]->edge() = new_inward_edges[i];
        inward_down_halfedges[i]->face() = new_faces[i == 0 ? degree - 1 : i - 1];

        inward_up_halfedges[i]->next() = orignal_halfedges[i];
        inward_up_halfedges[i]->twin() = inward_down_halfedges[i];
        inward_up_halfedges[i]->vertex() = new_vertices[i];
        inward_up_halfedges[i]->edge() = new_inward_edges[i];
        inward_up_halfedges[i]->face() = new_faces[i];

        face_cw_halfedges[i]->next() = inward_up_halfedges[i];
        face_cw_halfedges[i]->twin() = face_ccw_halfedges[i];
        face_cw_halfedges[i]->vertex() = new_vertices[i + 1 == degree ? 0 : i + 1];
        face_cw_halfedges[i]->edge() = new_face_edges[i];
        face_cw_halfedges[i]->face() = new_faces[i];

        face_ccw_halfedges[i]->next() = face_ccw_halfedges[i + 1 == degree ? 0 : i + 1];
        face_ccw_halfedges[i]->twin() = face_cw_halfedges[i];
        face_ccw_halfedges[i]->vertex() = new_vertices[i];
        face_ccw_halfedges[i]->edge() = new_face_edges[i];
        face_ccw_halfedges[i]->face() = f;

        orignal_halfedges[i]->next() = inward_down_halfedges[i + 1 == degree ? 0 : i + 1];
        orignal_halfedges[i]->face() = new_faces[i];

        new_inward_edges[i]->halfedge() = inward_down_halfedges[i];
        new_face_edges[i]->halfedge() = face_cw_halfedges[i];
        new_vertices[i]->halfedge() = face_ccw_halfedges[i];
        new_vertices[i]->pos = orignal_vertices[i]->pos;
        new_faces[i]->halfedge() = face_cw_halfedges[i];
    }
    f->halfedge() = face_ccw_halfedges[0];
    return f;
}

/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    // (void)new_halfedges;
    // (void)start_positions;
    // (void)face;
    // (void)tangent_offset;
    // (void)normal_offset;

    int number_of_vertices = (int)new_halfedges.size();
    Vec3 sum_of_all_vertices;
    for(int i = 0; i < number_of_vertices; i++) {
        sum_of_all_vertices = sum_of_all_vertices + start_positions[i];
    }
    Vec3 center = sum_of_all_vertices / number_of_vertices;

    for(int i = 0; i < (int)new_halfedges.size(); i++) {
        Vec3 face_normal_direction = face->normal();
        Vec3 normal_direction_delta = normal_offset * face_normal_direction;

        float percentage = 0.5;
        percentage = percentage + tangent_offset;

        if(percentage > 0.95) percentage = 0.95;
        if(percentage < 0.05) percentage = 0.05;
        Vec3 new_postion =
            center * (percentage) + start_positions[i] * (1 - percentage) + normal_direction_delta;
        new_halfedges[i]->vertex()->pos = new_postion;
    }
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {
    printf("triangulating\n");
    size_t number_of_faces = n_faces();
    Halfedge_Mesh::FaceRef f = faces_begin();
    // for(Halfedge_Mesh::FaceRef f = faces_begin(); f < faces_end(); f++) {

    while (number_of_faces)
    {
        number_of_faces--;
        printf("face degree is %d\n", f->degree());
        if(f->degree() == 3 || f->is_boundary()) {
            printf("face degree is 3 or boundary\n");
            f++;
            continue;
        }
        size_t n = f->degree() - 3;
        std::vector<Halfedge_Mesh::EdgeRef> edges;
        std::vector<Halfedge_Mesh::FaceRef> faces;
        std::vector<Halfedge_Mesh::HalfedgeRef> halfedges;
        
        // allocate memory for new edges, faces, and halfedges
        for(size_t i = 0; i < n; i++) {
            edges.push_back(new_edge());
            faces.push_back(new_face());
            halfedges.push_back(new_halfedge());
            halfedges.push_back(new_halfedge());
        }

        Halfedge_Mesh::HalfedgeRef halfedge0 = f->halfedge();
        
        // the first traingle
        // edges[0]->halfedge() = halfedges[0];
        
        // halfedge0->next()->next() = halfedges[0];
        halfedges[0]->Halfedge_Mesh::Halfedge::set_neighbors(halfedge0, halfedges[1], halfedge0->next()->next()->vertex(), edges[0], f);
        
        Halfedge_Mesh::HalfedgeRef h = halfedge0->next()->next();
        for(size_t i = 0; i < n - 1; i++)
        {
            halfedges[2 * i + 1]->Halfedge_Mesh::Halfedge::set_neighbors(h, halfedges[2 * i], halfedge0->vertex(), edges[i], faces[i]);
            h++;
            halfedges[2 * i + 2]->Halfedge_Mesh::Halfedge::set_neighbors(halfedges[2 * i + 1], halfedges[2 * i + 3], h->vertex(), edges[i + 1], faces[i]);
            faces[i]->halfedge() = halfedges[2 * i + 1];
            edges[i]->halfedge() = halfedges[2 * i + 1];
        }

        halfedges[2 * n - 1]->Halfedge_Mesh::Halfedge::set_neighbors(h, halfedges[2 * n - 2], halfedge0->vertex(), edges[n - 1], faces[n - 1]);
        faces[n - 1]->halfedge() = halfedges[2 * n - 1];
        edges[n - 1]->halfedge() = halfedges[2 * n - 1];
        h->face() = faces[n - 1];
        h->next()->face() = faces[n - 1];
        h->next()->next() = halfedges[2 * n - 1];

        h = halfedge0->next();
        Halfedge_Mesh::HalfedgeRef h_next = h->next();
        h->next() = halfedges[0];
        h = h_next;
        for(size_t i = 0; i < n - 1; i++)
        {
            h_next = h->next();
            h->face() = faces[i];
            h->next() = halfedges[2 * i + 2];
            h = h_next;
        }
        
        f++;
    }
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        v->new_pos = v->pos;
    }
    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        e->new_pos = e->center();
    }
    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        f->new_pos = f->center();
    }
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces
    for(FaceRef f = faces_begin(); f != faces_end(); f++) {
        f->new_pos = f->center();
    }
    // Edges
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Vec3 v1 = e->halfedge()->vertex()->pos;
        Vec3 v2 = e->halfedge()->twin()->vertex()->pos;
        Vec3 f1 = e->halfedge()->face()->new_pos;
        Vec3 f2 = e->halfedge()->twin()->face()->new_pos;
        e->new_pos = (v1 + v2 + f1 + f2) / 4;
    }
    // Vertices
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        Vec3 m, f, m_bar, f_bar;
        Vec3 p = v->pos;
        int n = 0;

        HalfedgeRef h0 = v->halfedge();
        HalfedgeRef h = h0;
        while(true) {
            m += h->edge()->new_pos;
            f += h->face()->new_pos;
            n++;
            h = h->twin()->next();
            if(h == h0) break;
        }
        m_bar = m / n;
        f_bar = f / n;
        v->new_pos = (f_bar / n) + (2 * m_bar / n) + (p * (n - 3) / n);
    }
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        //calculate n and u
        int n = v->degree();
        float u;
        if (n==3){
            u = 3.0/16.0;
        }
        else{
            u = 3.0/(8.0*n);
        }

        //calculate new position for all the vertices
        HalfedgeRef he0 = v->halfedge();
        HalfedgeRef he = he0;
        Vec3 new_postion = (1 - n*u) * he0->vertex()->pos;
        do{
            he = he->twin();
            new_postion += he -> vertex() -> pos * u;
            he = he->next();
        } while (he != he0);
        v->new_pos = new_postion;
    }
    // Next, compute the updated vertex positions associated with edges.
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Vec3 new_position;
        HalfedgeRef he = e->halfedge();
        new_position += (3.0 / 8.0) * he->vertex()->pos;
        new_position += (3.0 / 8.0) * he->twin()->vertex()->pos;
        new_position += (1.0 / 8.0) * he->next()->twin()->vertex()->pos;
        new_position += (1.0 / 8.0) * he->twin()->next()->twin()->vertex()->pos;
        e->new_pos = new_position;
    }
    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)
    int n = n_edges();
    EdgeRef e = edges_begin();
    for(int i = 0; i < n; i++) {

        // get the next edge NOW!
        EdgeRef nextEdge = e;
        nextEdge++;

        // now, even if splitting the edge deletes it...
        if(!e->is_new) {
            VertexRef new_v = split_edge(e).value();
            new_v -> new_pos = e->new_pos;
        }

        // ...we still have a valid reference to the next edge.
        e = nextEdge;
    }
    // Finally, flip any new edge that connects an old and new vertex.
    for(EdgeRef e = edges_begin(); e != edges_end(); e++) {
        if (e ->is_new){
            bool vertex1_new = e->halfedge()->vertex()->is_new;
            bool vertex2_new = e->halfedge()->twin()->vertex()->is_new;
            if (vertex1_new != vertex2_new){
                flip_edge(e);
            }
        }
    }
    // Copy the updated vertex positions to the subdivided mesh.
    for(VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        v->pos = v->new_pos;
    }
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.

        // solves for the optimal midpoint position

        Halfedge_Mesh::VertexRef vertex0 = e->halfedge()->vertex();
        Halfedge_Mesh::VertexRef vertex1 = e->halfedge()->twin()->vertex();
        Mat4 K = vertex_quadrics[vertex0] + vertex_quadrics[vertex1];
        
        float cost0 = dot(Vec4(vertex0->pos, 1.0), (K * Vec4(vertex0->pos, 1.0)));
        float cost1 = dot(Vec4(vertex1->pos, 1.0), (K * Vec4(vertex1->pos, 1.0)));
        if(cost0 < cost1) {
            optimal = vertex0->pos;
        } else {
            optimal = vertex1->pos;
        }
        cost = dot(Vec4(optimal, 1.0), (K * Vec4(optimal, 1.0)));
    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    for(Halfedge_Mesh::FaceRef f = faces_begin(); f != faces_end(); f++)
    {
        Vec3 N = f->normal();
        Vec3 p = f->center();
        float d = -dot(N, p);
        Vec4 u(p, 1);
        Vec4 v(N, d);
        face_quadrics[f] = outer(v, v);
    }

    for(Halfedge_Mesh::VertexRef v = vertices_begin(); v != vertices_end(); v++) {
        Halfedge_Mesh::HalfedgeRef halfedge0 = v->halfedge();
        Halfedge_Mesh::HalfedgeRef h = halfedge0->twin()->next();
        
        Mat4 face_quadric_sum = face_quadrics[halfedge0->face()];
        while (h != halfedge0) {
            face_quadric_sum += face_quadrics[h->face()];
            h = h->twin()->next();
        }
        vertex_quadrics[v] = face_quadric_sum;
    }
    
    for(Halfedge_Mesh::EdgeRef e = edges_begin(); e != edges_end(); e++) {
        Edge_Record er(vertex_quadrics, e);
        edge_records[e] = er;
        edge_queue.insert(er);
    }

    size_t target_face_number = (size_t)(faces.size() / 4);
    if(target_face_number <= 0) return false;


    while(target_face_number-- && !edge_queue.queue.empty())
    {
        Edge_Record bestEdge = edge_queue.top();
        edge_queue.pop();
        
        Mat4 q1 = vertex_quadrics[bestEdge.edge->halfedge()->vertex()];
        Mat4 q2 = vertex_quadrics[bestEdge.edge->halfedge()->twin()->vertex()];
        Mat4 q_sum = q1 + q2;
        vertex_quadrics.erase(bestEdge.edge->halfedge()->vertex());
        vertex_quadrics.erase(bestEdge.edge->halfedge()->twin()->vertex());


        std::unordered_set<EdgeRef> neighbors;
        Halfedge_Mesh::HalfedgeRef halfedge0 = bestEdge.edge->halfedge();
        Halfedge_Mesh::HalfedgeRef h = halfedge0->twin()->next();
        // neighbors.insert(halfedge0->edge());
        while (h != halfedge0) {
            neighbors.insert(h->edge());
            h = h->twin()->next();
        }
        Halfedge_Mesh::HalfedgeRef halfedge1 = bestEdge.edge->halfedge()->twin();
        h = halfedge1->twin()->next();
        // neighbors.insert(halfedge1->edge());
        while (h != halfedge1) {
            neighbors.insert(h->edge());
            h = h->twin()->next();
        }
        
        std::unordered_set<EdgeRef>::iterator it;
        for (it = neighbors.begin(); it != neighbors.end(); it++) {
            if (edge_records.find(*it) != edge_records.end()) {
                edge_queue.remove(edge_records[*it]);
            }
        }
        
        auto optional_vertex = collapse_edge_erase(bestEdge.edge);
        edge_records.erase(bestEdge.edge);
        std::unordered_set<EdgeRef> new_neighbors;
        if (optional_vertex.has_value()) {
            Halfedge_Mesh::VertexRef v = optional_vertex.value();
            v->pos = bestEdge.optimal; // the optimal point, if we were to collapse this edge next
            vertex_quadrics[v] = q_sum;
            Halfedge_Mesh::HalfedgeRef halfedge0 = v->halfedge();
            Halfedge_Mesh::HalfedgeRef h = halfedge0->twin()->next();
            new_neighbors.insert(halfedge0->edge());
            while (h != halfedge0) {
                new_neighbors.insert(h->edge());
                h = h->twin()->next();
            }
        }
        
        if (optional_vertex.has_value()) {
            for (it = new_neighbors.begin(); it != new_neighbors.end(); it++) {
                if (edge_records.find(*it) != edge_records.end()) {
                    Edge_Record er(vertex_quadrics, *it);
                    edge_records[*it] = er;
                    edge_queue.insert(er);
                }
            }
        }
        else
        {
            for (it = neighbors.begin(); it != neighbors.end(); it++) {
                if (edge_records.find(*it) != edge_records.end()) {
                    Edge_Record er(vertex_quadrics, *it);
                    edge_records[*it] = er;
                    edge_queue.insert(er);
                }
            }
        }
    }

    return true;
    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.
}