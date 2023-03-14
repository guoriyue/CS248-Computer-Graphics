#define WIDTH 4
// .data()
// export void build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
export void SAH(uniform size_t idx, uniform size_t max_leaf_size) {
    // nodes.clear();
    // primitives = std::move(prims);

    // BBox bb;
    // for(size_t i = 0; i < primitives.size(); ++i) {
    //     bb.enclose(primitives[i].bbox());
    // }
    // // return 0, should be the index of the root node
    // size_t root_node_addr = new_node();
    // Node& node = nodes[root_node_addr];
    // node.bbox = bb;
    // node.start = 0;
    // node.size = primitives.size();
    // SAH(root_node_addr, max_leaf_size);
}

// task void SAH(size_t idx, size_t max_leaf_size) {
//     if (nodes[idx].size <= max_leaf_size) {
//         return;
//     }
//     // Create bounding boxes for children
//     BBox split_leftBox;
//     BBox split_rightBox;
//     BBox bbox = nodes[idx].bbox;
//     int bucket_num = 8;

//     int min_cost_axis = 0;
//     float min_cost = 0x7fffffff;
//     float min_cost_split = 0;


//     size_t rangel = 0;
//     size_t ranger = 0;
//     // can't use buckets[3][bucket_num], will cause bus fault
//     for (int i = 0; i < 3; i++) {
//         BBox buckets[bucket_num];
//         int prim_count[bucket_num];
//         for (int i = 0; i < 3; i++) {
//             for (int j = 0; j < bucket_num; j++) {
//                 prim_count[j] = 0;
//             }
//         }
//         float my_min = bbox.min[i];
//         float my_max = bbox.max[i];
//         float interval = (my_max - my_min) / (bucket_num + 0.0);
//         for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
//             Primitive& p = primitives[j];
//             BBox pbb = p.bbox();
            
//             int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
//             bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);
//             buckets[bucket_idx].enclose(pbb);
//             prim_count[bucket_idx]++;
//         }
        
        
//         for (int j = 0; j < bucket_num - 1; j++) {
//             float left_surface_area = 0;
//             float right_surface_area = 0;
//             float left_prim_count = 0;
//             float right_prim_count = 0;

//             BBox left_bbox;
//             BBox right_bbox;
    
//             for (int k = 0; k <= j; k++) {
//                 left_bbox.enclose(buckets[k]);
//                 left_prim_count += prim_count[k];
//             }
//             for (int k = j + 1; k < bucket_num; k++) {
//                 right_bbox.enclose(buckets[k]);
//                 right_prim_count += prim_count[k];
//             }

//             left_surface_area = left_bbox.surface_area();
//             right_surface_area = right_bbox.surface_area();

//             float total_cost = left_surface_area * left_prim_count + right_surface_area * right_prim_count;
//             if (total_cost < min_cost) {
//                 min_cost = total_cost;
//                 min_cost_axis = i;
//                 min_cost_split = my_min + (j + 1) * interval;
//                 split_leftBox = left_bbox;
//                 split_rightBox = right_bbox;
//                 rangel = left_prim_count;
//                 ranger = right_prim_count;
//             }
//         }
//     }

//     // need to reorganize primitives so that the children are contiguous ranges of primitives
//     int first = nodes[idx].start;
//     for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) 
//     {
//         Primitive& p = primitives[j];
//         BBox pbb = p.bbox();
//         if (pbb.center()[min_cost_axis] < min_cost_split)
//         {
//             std::swap(primitives[j], primitives[first]);
//             ++first;
//         }
//     }

//     size_t startl = nodes[idx].start;
//     size_t startr = startl + rangel;

//     // create child nodes
//     size_t node_addr_l = new_node();
//     size_t node_addr_r = new_node();
//     nodes[idx].l = node_addr_l;
//     nodes[idx].r = node_addr_r;

//     nodes[node_addr_l].bbox = split_leftBox;
//     nodes[node_addr_l].start = startl;
//     nodes[node_addr_l].size = rangel;

//     nodes[node_addr_r].bbox = split_rightBox;
//     nodes[node_addr_r].start = startr;
//     nodes[node_addr_r].size = ranger;

//     SAH(node_addr_l, max_leaf_size);
//     SAH(node_addr_r, max_leaf_size);
// }