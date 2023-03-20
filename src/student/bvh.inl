
#include "../rays/bvh.h"

// #include "../rays/Kmeans.h"
#include <stack>
// #include "../lib/mathlib.h"
#include "../rays/ispc_bvh.h"
#define SIMD_WIDTH 4
// using namespace std;
  //聚类数
static const size_t K=SIMD_WIDTH;  

// struct KmeansBuildData {
//   KmeansBuildData(BBox bb,KBVHNode **dst)
//       : bb(bb), node(dst) {}
//   BBox bb;         ///  包围片元的包围盒
//   KBVHNode **node;  ///  对应的BVHNODE
// };


namespace PT {

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    printf("Inside BVH::build\n");

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Keep these two lines of code in your solution. They clear the list of nodes and
    // initialize member variable 'primitives' as a vector of the scene prims
    nodes.clear();
    primitives = std::move(prims);

    // // TODO (PathTracer): Task 3
    // // Modify the code ahead to construct a BVH from the given vector of primitives and maximum leaf
    // // size configuration.
    // //
    // // Please use the SAH as described in class.  We recomment the binned build from lecture.
    // // In general, here is a rough sketch:
    // //
    // //  For each axis X,Y,Z:
    // //     Try possible splits along axis, evaluate SAH for each
    // //  Take minimum cost across all axes.
    // //  Partition primitives into a left and right child group
    // //  Compute left and right child bboxes
    // //  Make the left and right child nodes.
    // //
    // //
    // // While a BVH is conceptually a tree structure, the BVH class uses a single vector (nodes)
    // // to store all the nodes. Therefore, BVH nodes don't contain pointers to child nodes,
    // // but rather the indices of the
    // // child nodes in this array. Hence, to get the child of a node, you have to
    // // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // // returns the index of a newly added node.
    // //
    // // As an example of how to make nodes, the starter code below builds a BVH with a
    // // root node that encloses all the primitives and its two descendants at Level 2.
    // // For now, the split is hardcoded such that the first primitive is put in the left
    // // child of the root, and all the other primitives are in the right child.
    // // There are no further descendants.

    // // edge case
    // if(primitives.empty()) {
    //     return;
    // }

    // // compute bounding box for all primitives
    // BBox bb;
    // for(size_t i = 0; i < primitives.size(); ++i) {
    //     bb.enclose(primitives[i].bbox());
    // }

    // // set up root node (root BVH). Notice that it contains all primitives.
    // size_t root_node_addr = new_node();
    // Node& node = nodes[root_node_addr];
    // node.bbox = bb;
    // node.start = 0;
    // node.size = primitives.size();

    // // Create bounding boxes for children
    // BBox split_leftBox;
    // BBox split_rightBox;

    // // compute bbox for left child
    // Primitive& p = primitives[0];
    // BBox pbb = p.bbox();
    // split_leftBox.enclose(pbb);

    // // compute bbox for right child
    // for(size_t i = 1; i < primitives.size(); ++i) {
    //     Primitive& p = primitives[i];
    //     BBox pbb = p.bbox();
    //     split_rightBox.enclose(pbb);
    // }

    // // Note that by construction in this simple example, the primitives are
    // // contiguous as required. But in the students real code, students are
    // // responsible for reorganizing the primitives in the primitives array so that
    // // after a SAH split is computed, the chidren refer to contiguous ranges of primitives.

    // size_t startl = 0;  // starting prim index of left child
    // size_t rangel = 1;  // number of prims in left child
    // size_t startr = startl + rangel;  // starting prim index of right child
    // size_t ranger = primitives.size() - rangel; // number of prims in right child

    // // create child nodes
    // size_t node_addr_l = new_node();
    // size_t node_addr_r = new_node();
    // nodes[root_node_addr].l = node_addr_l;
    // nodes[root_node_addr].r = node_addr_r;

    // nodes[node_addr_l].bbox = split_leftBox;
    // nodes[node_addr_l].start = startl;
    // nodes[node_addr_l].size = rangel;

    // nodes[node_addr_r].bbox = split_rightBox;
    // nodes[node_addr_r].start = startr;
    // nodes[node_addr_r].size = ranger;


    BBox bb;
    for(size_t i = 0; i < primitives.size(); ++i) {
        bb.enclose(primitives[i].bbox());
    }
    // return 0, should be the index of the root node
    size_t root_node_addr = new_node4();
    Node4& node = nodes[root_node_addr];
    // size_t root_node_addr = new_node8();
    // Node8& node = nodes[root_node_addr];
    // size_t root_node_addr = new_node();
    // Node& node = nodes[root_node_addr];
    
    node.bbox = bb;
    node.start = 0;
    node.size = primitives.size();
    
    // SAH(root_node_addr, max_leaf_size);

    // SAH4(root_node_addr, max_leaf_size);






//     //构造函数  
// template<typename Primitive>
// void Kmeans<Primitive>::InitKmeans(size_t iterCount, size_t K, size_t P,std::vector<Primitive> primitives){
//    kmeans = new kmeans(m_iterations,m_K,m_P,primitives);
printf("before init kmeans\n");
kmeans_root = new Kmeans(2,SIMD_WIDTH,5,std::vector<size_t>(0, primitives.size()));
printf("before 11232 init kmeans\n");
   kmeans_root->m_iterations=2;
   printf("befor444 e init kmeans\n");
   kmeans_root->m_K=SIMD_WIDTH;
   printf("before 555 init kmeans\n");
   kmeans_root->m_P=5;
//    kmeans.primitives=primitives;
    // kmeans.indexOfPrimitives=new size_t[primitives.size()];
    for(size_t i=0;i<primitives.size();++i){
      kmeans_root->indexOfPrimitives.push_back(i);
    }
   kmeans_root->cluster=new Cluster[kmeans_root->m_K];
   kmeans_root->children=new Kmeans* [kmeans_root->m_K];
   printf("after init kmeans\n");
   BBox world;
   for(size_t i=0;i<primitives.size();++i){
      world.enclose(primitives[i].bbox());
    } 
   
   /*
   //初始点随机选取
   vector <Vector3D> kCentroids=getRandCentroids(world.min,world.max,m_K,5);
    for(size_t i=0;i<m_K;++i){
      cluster[i].representive=BBox(kCentroids[i]);
  }*/
  
  auto length=[&](Vec3 v){
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  };
  auto getRandCentroidsOnMesh=[&](size_t k,size_t p,std::vector<size_t> indexOfPrimitives){
    // std::vector <BBox> centroids;
    // for(size_t i=0;i<k;++i){
    //   BBox centroid;
    //   for(size_t j=0;j<p;++j){
    //     size_t index=rand()%primitives.size();
    //     centroid.enclose(primitives[index].bbox());
    //   }
    //   centroids.push_back(centroid);
    // }
    // return centroids;

   std::vector <BBox> kCentroids;
   size_t idx_primitive;
   
   //第一个随机点
//    srand((unsigned)time(NULL));
    size_t prim_size = indexOfPrimitives[indexOfPrimitives.size()-1] - indexOfPrimitives[0] + 1;
   idx_primitive=rand()%prim_size + indexOfPrimitives[0];
   kCentroids.push_back(primitives[idx_primitive].bbox());
   
   printf("in 111 getRandCentroidsOnMesh\n");
    // srand((unsigned)time(NULL));
   //选取之后k-1个点
   for(unsigned long i=1;i<k;++i){
     //随机p个点
      std::vector<BBox> tempP;
      for(unsigned long j=0;j<p;j++){
        // idx_primitive=rand()%primitives.size();
        // tempP.push_back(primitives[idx_primitive].bbox());
        idx_primitive=rand()%prim_size + indexOfPrimitives[0];
        tempP.push_back(primitives[idx_primitive].bbox());	
      }
    
      printf("in 222 getRandCentroidsOnMesh\n");
      //选取与之前算出的点距离最远的点作为下一个representive
      int index=0;
      double maxDistance=-1.0f;
      
      for(unsigned long k=0;k<p;k++)
      {
	for(unsigned long q=0;q<kCentroids.size();q++){
	   BBox bb;
	   bb.enclose(tempP[k]);
	   bb.enclose(kCentroids[q]);
	    double distance=length(bb.max - bb.min);
	    if(distance>maxDistance){
	       maxDistance=distance;
	       index=k;
	    }
	}	
      }
     printf("in 333 getRandCentroidsOnMesh\n");
      kCentroids.push_back(tempP[index]);
  }
  return kCentroids;

  };
  printf("before getRandCentroidsOnMesh\n");
  
   std::vector <BBox> kCentroids=getRandCentroidsOnMesh(kmeans_root->m_K,kmeans_root->m_P, kmeans_root->indexOfPrimitives);
   printf("after getRandCentroidsOnMesh\n");
   for(size_t i=0;i<kmeans_root->m_K;++i){
     kmeans_root->cluster[i].representive=kCentroids[i];
     printf("in initialization cluster[i].indexOfPrimitives.size()=%lu\n",kmeans_root->cluster[i].indexOfPrimitives.size());
     printf("kmeans->cluster[i].representive.min.x=%f\n",kmeans_root->cluster[i].representive.min.x);
        printf("kmeans->cluster[i].representive.min.y=%f\n",kmeans_root->cluster[i].representive.min.y);
        printf("kmeans->cluster[i].representive.min.z=%f\n",kmeans_root->cluster[i].representive.min.z);
        printf("kmeans->cluster[i].representive.max.x=%f\n",kmeans_root->cluster[i].representive.max.x);
        printf("kmeans->cluster[i].representive.max.y=%f\n",kmeans_root->cluster[i].representive.max.y);
        printf("kmeans->cluster[i].representive.max.z=%f\n",kmeans_root->cluster[i].representive.max.z);
  }
printf("max_leaf_size=%lu\n",max_leaf_size);
constructKaryTree(kmeans_root,max_leaf_size);
printf("finish BVH::build\n");
// print();

}

template<typename Primitive>
void BVH<Primitive>::allocateKaryTree(Kmeans* kmeans,size_t maxLeafNum)
{
    printf("start allocateKaryTree::build\n");
    auto length=[&](Vec3 v){
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  };
auto calDistance=[&](BBox b1, BBox b2){
  
 
    double min_value=length(b1.min-b2.min);
    double max_value=length(b1.max-b2.max);
    double res=min_value+max_value;
    return res;
    
    
 //  return (b1.centroid()-b2.centroid()).norm2();
};


// calDistance(BBox b1, BBox b2){
  
 
//     double min_value=length(b1.min-b2.min);
//     double max_value=length(b1.max-b2.max);
//     double res=min_value+max_value;
//     return res;
    
    
//  //  return (b1.centroid()-b2.centroid()).norm2();
// }


  for(size_t iter=0;iter<kmeans->m_iterations;++iter){
    //更新representives cluster重新分配
    for(size_t i=0;i<kmeans->m_K;++i){
      
      if(iter!=0){
       kmeans->cluster[i].updateRepresentive(); 
      }
      kmeans->cluster[i].reset();
      kmeans->cluster[i].indexOfPrimitives.clear();
    }
        
    //计算距离 分配至最近的cluster
    for(size_t idx_primitives=0;idx_primitives<kmeans->indexOfPrimitives.size();++idx_primitives){
         //用以记录最近的cluster
         size_t index=0;
	 double minDistance=std::numeric_limits< double >::max();
	 BBox temp=primitives[kmeans->indexOfPrimitives[idx_primitives]].bbox();
     printf("idx_primitives %lu\n",idx_primitives);
     printf("temp %lf %lf %lf\n",temp.min.x,temp.min.y,temp.min.z);
        printf("temp %lf %lf %lf\n",temp.max.x,temp.max.y,temp.max.z);

	 for(size_t idx_clusters=0;idx_clusters<kmeans->m_K;++idx_clusters){
        printf("kmeans->cluster[idx_clusters].representive %lf %lf %lf\n",kmeans->cluster[idx_clusters].representive.min.x,kmeans->cluster[idx_clusters].representive.min.y,kmeans->cluster[idx_clusters].representive.min.z);
        printf("kmeans->cluster[idx_clusters].representive %lf %lf %lf\n",kmeans->cluster[idx_clusters].representive.max.x,kmeans->cluster[idx_clusters].representive.max.y,kmeans->cluster[idx_clusters].representive.max.z);    	    
	     double dist=calDistance(temp,kmeans->cluster[idx_clusters].representive);
         printf("dist %lf\n",dist);
         printf("index %lu\n",index);
	     if(dist<minDistance){
	       minDistance=dist;
	       index=idx_clusters;
	    }
	}
	//allocate to correct cluster;
	kmeans->cluster[index].add(idx_primitives,primitives[kmeans->indexOfPrimitives[idx_primitives]].bbox());
    }    
  }
  printf("end allocateKaryTree::build\n");
}

template<typename Primitive>
void BVH<Primitive>::constructKaryTree(Kmeans* kmeans,size_t maxLeafNum)
{
  //构造本层结构
  printf("start constructKaryTree::build\n");
  allocateKaryTree(kmeans,maxLeafNum);
  //对本层中的每个cluster循环构造下层
  for(size_t i=0;i<kmeans->m_K;i++){
     //叶子cluster 该cluster[i]对应的children为NULL
     if(kmeans->cluster[i].indexOfPrimitives.size()<maxLeafNum*kmeans->m_K){
       kmeans->children[i]=NULL;
       continue;
     }
     if(kmeans->cluster[i].representive.min.x>0x3fffff){
        kmeans->children[i]=NULL;
        continue;
     }
     
     //否则DFS
     std::vector <size_t> pTemp;
     printf("cluster[i].indexOfPrimitives.size()=%lu\n",kmeans->cluster[i].indexOfPrimitives.size());
     for(size_t p=0;p<kmeans->cluster[i].indexOfPrimitives.size();++p){
          pTemp.push_back(kmeans->cluster[i].indexOfPrimitives[p]);
     }
     kmeans->children[i]=new Kmeans(kmeans->m_iterations,kmeans->m_K,kmeans->m_P,pTemp);
     
     constructKaryTree(kmeans->children[i], maxLeafNum);
  }     
  printf("end constructKaryTree::build\n");
}
// template<typename Primitive>
// Kmeans* BVH<Primitive>::new_Kmeans(size_t iterations,size_t K,size_t P,std::vector <size_t> indexOfPrimitives)
// {
//   Kmeans* kmeans=new Kmeans(iterations,K,P,indexOfPrimitives);
//   return kmeans;
// }


// Kmeans<Primitive>* k=new Kmeans(2,8,5,primitives);
// // timer.start();
// k->constructKaryTree();
// //k->run();
// // timer.stop();
// // cout<<"K-aryTree 自顶向下构建共耗时(sec):  "<<timer.duration()<<endl;

// // timer.start();
// k->buttom2Top();
// KBVHNode* root=k->root;
// timer.stop();
// cout<<"自底向上构建共耗时(sec):  "<<timer.duration()<<endl;

// SAH8(root_node_addr, max_leaf_size);
// printf("SAH8 done\n");
// }

template<typename Primitive>
void BVH<Primitive>::SAH(const size_t & idx, const size_t & max_leaf_size) {
    if (nodes[idx].size <= max_leaf_size) {
        return;
    }
    // Create bounding boxes for children
    BBox split_leftBox;
    BBox split_rightBox;
    BBox bbox = nodes[idx].bbox;
    int bucket_num = 8;

    int min_cost_axis = 0;
    float min_cost = 0x7fffffff;
    float min_cost_split = 0;


    size_t rangel = 0;
    size_t ranger = 0;
    // can't use buckets[3][bucket_num], will cause bus fault
    for (int i = 0; i < 3; i++) {
        BBox buckets[bucket_num];
        int prim_count[bucket_num];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < bucket_num; j++) {
                prim_count[j] = 0;
            }
        }
        float my_min = bbox.min[i];
        float my_max = bbox.max[i];
        float interval = (my_max - my_min) / (bucket_num + 0.0);
        for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
            Primitive& p = primitives[j];
            BBox pbb = p.bbox();
            
            int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
            bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);
            
            buckets[bucket_idx].enclose(pbb);
            prim_count[bucket_idx]++;
        }
        
        
        for (int j = 0; j < bucket_num - 1; j++) {
            float left_surface_area = 0;
            float right_surface_area = 0;
            float left_prim_count = 0;
            float right_prim_count = 0;

            BBox left_bbox;
            BBox right_bbox;
    
            for (int k = 0; k <= j; k++) {
                left_bbox.enclose(buckets[k]);
                left_prim_count += prim_count[k];
            }
            for (int k = j + 1; k < bucket_num; k++) {
                right_bbox.enclose(buckets[k]);
                right_prim_count += prim_count[k];
            }

            left_surface_area = left_bbox.surface_area();
            right_surface_area = right_bbox.surface_area();

            float total_cost = left_surface_area * left_prim_count + right_surface_area * right_prim_count;
            if (total_cost < min_cost) {
                min_cost = total_cost;
                min_cost_axis = i;
                min_cost_split = my_min + (j + 1) * interval;
                split_leftBox = left_bbox;
                split_rightBox = right_bbox;
                rangel = left_prim_count;
                ranger = right_prim_count;
            }
        }
    }

    // need to reorganize primitives so that the children are contiguous ranges of primitives
    int first = nodes[idx].start;
    for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) 
    {
        Primitive& p = primitives[j];
        BBox pbb = p.bbox();
        if (pbb.center()[min_cost_axis] < min_cost_split)
        {
            std::swap(primitives[j], primitives[first]);
            ++first;
        }
    }

    size_t startl = nodes[idx].start;
    size_t startr = startl + rangel;

    // create child nodes
    size_t node_addr_l = new_node();
    size_t node_addr_r = new_node();
    nodes[idx].l = node_addr_l;
    nodes[idx].r = node_addr_r;

    nodes[node_addr_l].bbox = split_leftBox;
    nodes[node_addr_l].start = startl;
    nodes[node_addr_l].size = rangel;

    nodes[node_addr_r].bbox = split_rightBox;
    nodes[node_addr_r].start = startr;
    nodes[node_addr_r].size = ranger;

    printf("node %lu, l %lu, r %lu, size %lu\n", idx, node_addr_l, node_addr_r, nodes[idx].size);
    printf("l %lu, r %lu\n", rangel, ranger);
    printf("bbox %f %f %f %f %f %f\n", bbox.min[0], bbox.min[1], bbox.min[2], bbox.max[0], bbox.max[1], bbox.max[2]);
    printf("l bbox %f %f %f %f %f %f\n", split_leftBox.min[0], split_leftBox.min[1], split_leftBox.min[2], split_leftBox.max[0], split_leftBox.max[1], split_leftBox.max[2]);
    printf("r bbox %f %f %f %f %f %f\n", split_rightBox.min[0], split_rightBox.min[1], split_rightBox.min[2], split_rightBox.max[0], split_rightBox.max[1], split_rightBox.max[2]);

    SAH(node_addr_l, max_leaf_size);
    SAH(node_addr_r, max_leaf_size);
}













template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    // Trace ret;
    // for(const Primitive& prim : primitives) {
    //     Trace hit = prim.hit(ray);
    //     ret = Trace::min(ret, hit);
    // }
    // return ret;
    // return find_hit(ray, 0);


    // Trace ret;
    // std::stack<size_t> node_stack;
    // node_stack.push(root_idx);

    // while(!node_stack.empty()) {
    //     auto node_idx = node_stack.top();
    //     node_stack.pop();
    //     auto& node = nodes[node_idx];

    //     Vec2 times0{};
    //     auto hit0 = node.bbox.hit(ray, times0);
    //     if(!hit0 || (ret.hit && ret.distance <= times0.x)) continue;

    //     if(node.is_leaf()) {
    //         auto node_end = node.start + node.size;
    //         for(auto i = node.start; i < node_end; ++i) {
    //             Trace hit = primitives[i].hit(ray);
    //             ret = Trace::min(ret, hit);
    //         }
    //     } else {
    //         Vec2 times1, times2;
    //         auto hit1 = nodes[node.l].bbox.hit(ray, times1);
    //         auto hit2 = nodes[node.r].bbox.hit(ray, times2);

    //         if(hit1 && hit2) {
    //             auto first = times1.x < times2.x ? node.l : node.r;
    //             auto second = times1.x < times2.x ? node.r : node.l;
    //             node_stack.push(first);
    //             node_stack.push(second);
    //         } else if(hit1) {
    //             node_stack.push(node.l);
    //         } else if(hit2) {
    //             node_stack.push(node.r);
    //         }
    //     }
    // }
    // return ret;

    return hit_kmeans(ray, kmeans_root);
    // return hit_queue4(ray);
    // return hit_queue8(ray);

    // queue for traversal
    // auto ispc_Vec3 = [](const Vec3& v) {
	// 	ispc::Vec3 res;
	// 	res.x = v.x; res.y = v.y; res.z = v.z;
	// 	return res;
	// };

    // Trace ret;
    // std::stack<size_t> node_stack;
    // node_stack.push(0);
    // while(!node_stack.empty()) {
    //     size_t node_idx = node_stack.top();
    //     node_stack.pop();
    //     const Node& node = nodes[node_idx];

    //     // with early return
    //     Vec2 times0{};
    //     bool hit0 = node.bbox.hit(ray, times0);
    //     if(!hit0 || (ret.hit && ret.distance <= times0.x))
    //     {
    //         continue;
    //     }
        
    //     if(node.is_leaf()) {
    //         size_t node_end = node.start + node.size;
    //         for(size_t i = node.start; i < node_end; ++i) {
    //             // printf("hit leaf node\n");
    //             Trace hit = primitives[i].hit(ray);
    //             ret = Trace::min(ret, hit);
    //         }
    //     } else {
    //         // std::vector<ispc::Vec2> ispc_times;
    //         // std::vector<ispc::Node> ispc_nodes;
    //         // ispc_times.resize(2);
    //         // ispc_nodes.resize(2);
    //         ispc::Vec2 *ispc_times = new ispc::Vec2[2];
    //         // ispc::Node *ispc_nodes = new ispc::Node[2];
            
    //         // ispc_nodes[0].bbox.min = ispc_Vec3(nodes[node.l].bbox.min);
    //         // ispc_nodes[0].bbox.max = ispc_Vec3(nodes[node.l].bbox.max);
    //         // ispc_nodes[1].bbox.min = ispc_Vec3(nodes[node.r].bbox.min);
    //         // ispc_nodes[1].bbox.max = ispc_Vec3(nodes[node.r].bbox.max);
    //         // ispc_nodes[0].start = nodes[node.l].start;
    //         // ispc_nodes[0].size = nodes[node.l].size;
    //         // ispc_nodes[0].l = nodes[node.l].l;
    //         // ispc_nodes[0].r = nodes[node.l].r;
    //         // ispc_nodes[1].start = nodes[node.r].start;
    //         // ispc_nodes[1].size = nodes[node.r].size;
    //         // ispc_nodes[1].l = nodes[node.r].l;
    //         // ispc_nodes[1].r = nodes[node.r].r;
            
    //         // ispc_nodes[0] = ispc::Node(ispc_Vec3(nodes[node.l].bbox.min), ispc_Vec3(nodes[node.l].bbox.max), nodes[node.l].start, nodes[node.l].size, nodes[node.l].l, nodes[node.l].r);
    //         // ispc_nodes[1] = ispc::Node(ispc_Vec3(nodes[node.r].bbox.min), ispc_Vec3(nodes[node.r].bbox.max), nodes[node.r].start, nodes[node.r].size, nodes[node.r].l, nodes[node.r].r);

    //         // ispc::Ray ispc_ray = ispc::Ray(ispc_Vec3(ray.point), ispc_Vec3(ray.dir));
    //         ispc::Ray ispc_ray;
    //         ispc_ray.point = ispc_Vec3(ray.point);
    //         ispc_ray.dir = ispc_Vec3(ray.dir);

    //         ispc::BBox *ispc_bbox = new ispc::BBox[2];
    //         ispc_bbox[0].min = ispc_Vec3(nodes[node.l].bbox.min);
    //         ispc_bbox[0].max = ispc_Vec3(nodes[node.l].bbox.max);
    //         ispc_bbox[1].min = ispc_Vec3(nodes[node.r].bbox.min);
    //         ispc_bbox[1].max = ispc_Vec3(nodes[node.r].bbox.max);
            
            
    //         bool* ispc_hits = new bool[2];
    //         ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);

    //         Vec2 times1, times2;
    //         times1.x = ispc_times[0].x;
    //         times1.y = ispc_times[0].y;
    //         times2.x = ispc_times[1].x;
    //         times2.y = ispc_times[1].y;
    //         // bool hit1 = nodes[node.l].bbox.hit(ray, times1);
    //         // bool hit2 = nodes[node.r].bbox.hit(ray, times2);
    //         bool hit1 = ispc_hits[0];
    //         bool hit2 = ispc_hits[1];
    //         if(hit1 && hit2) {
    //             size_t first = times1.x < times2.x ? node.l : node.r;
    //             size_t second = times1.x < times2.x ? node.r : node.l;
    //             node_stack.push(first);
    //             node_stack.push(second);
    //         } else if(hit1) {
    //             node_stack.push(node.l);
    //         } else if(hit2) {
    //             node_stack.push(node.r);
    //         }
    //     }
    // }
    // return ret;

    // ispc::Ray ispc_ray;
    // ispc_ray.point = ispc_Vec3(ray.point);
    // ispc_ray.dir = ispc_Vec3(ray.dir);
    // ispc::Node* ispc_nodes = new ispc::Node[nodes.size()];
    // for(int i = 0; i < nodes.size(); i++) {
    //     ispc_nodes[i].bbox.min = ispc_Vec3(nodes[i].bbox.min);
    //     ispc_nodes[i].bbox.max = ispc_Vec3(nodes[i].bbox.max);
    //     ispc_nodes[i].start = nodes[i].start;
    //     ispc_nodes[i].size = nodes[i].size;
    //     ispc_nodes[i].l = nodes[i].l;
    //     ispc_nodes[i].r = nodes[i].r;
    // }
    // ispc::Trace ispc_ret;

    // ispc::find_hit(ispc_ray, 0, ispc_nodes, ispc_ret);

    // Trace ret;
    // ret.hit = ispc_ret.hit;
    // ret.distance = ispc_ret.distance;
    // ret.position = Vec3(ispc_ret.position.x, ispc_ret.position.y, ispc_ret.position.z);
    // ret.normal = Vec3(ispc_ret.normal.x, ispc_ret.normal.y, ispc_ret.normal.z);
    // ret.origin = Vec3(ispc_ret.origin.x, ispc_ret.origin.y, ispc_ret.origin.z);
    // ret.material = ispc_ret.material;
    // return ret;
}



template<typename Primitive>
Trace BVH<Primitive>::hit_kmeans(const Ray& ray, Kmeans *kmeans) const {
    // Trace ret;
    // bool visit = false;
    // for(int j = 0; j < SIMD_WIDTH; j++) {
    //     if(kmeans->children[j] != NULL && kmeans->cluster[j].indexOfPrimitives.size()) {
    //         visit = true;

    //         Vec2 times;
    //         bool hit = kmeans->cluster[j].representive.hit(ray, times);
    //         if(hit) {
    //             Trace hit = hit_kmeans(ray, kmeans->children[j]);
    //             ret = Trace::min(ret, hit);
    //         }
    //     }
    // }
    // if(!visit)
    // {
    //     for(size_t p=0;p<kmeans->indexOfPrimitives.size();++p){
    //         Trace hit = primitives[kmeans->indexOfPrimitives[p]].hit(ray);
    //         ret = Trace::min(ret, hit);
    //     }
    // }

    // return ret;




    auto ispc_Vec3 = [](const Vec3& v) {
		ispc::Vec3 res;
		res.x = v.x; res.y = v.y; res.z = v.z;
		return res;
	};

    // auto ispc_Vec2 = [](const Vec2& v) {
    //     ispc::Vec2 res;
    //     res.x = v.x; res.y = v.y;
    //     return res;
    // };

    // ispc::Ray ispc_ray;
    // ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    // ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];
    // bool* ispc_hits = new bool[SIMD_WIDTH];
    // Vec2 *times = new Vec2[SIMD_WIDTH];
    // ispc::Trace *ispc_ret = new ispc::Trace[1];
    // ispc::Triangle* ispc_triangles = new ispc::Triangle[primitives.size()];
    // ispc::Node* ispc_nodes = new ispc::Node[nodes.size()];
    // for(int i = 0; i < primitives.size(); i++) {
    //     // ispc_triangles[i].v0 = ispc_Vec3(primitives[i].v0);
    //     // ispc_triangles[i].v1 = ispc_Vec3(primitives[i].v1);
    //     // ispc_triangles[i].v2 = ispc_Vec3(primitives[i].v2);
    //     // ispc_triangles[i].normal = ispc_Vec3(primitives[i].normal);
    //     // ispc_triangles[i].material = primitives[i].material;
    //     ispc_triangles[i].v_0 = vertex_list[primitives[i].v0];
    //     ispc_triangles[i].v_1 = vertex_list[primitives[i].v1];
    //     ispc_triangles[i].v_2 = vertex_list[primitives[i].v2];
    // }
    // for(int i = 0; i < nodes.size(); i++) {
    //     ispc_nodes[i].bbox.min = ispc_Vec3(nodes[i].bbox.min);
    //     ispc_nodes[i].bbox.max = ispc_Vec3(nodes[i].bbox.max);
    //     ispc_nodes[i].start = nodes[i].start;
    //     ispc_nodes[i].size = nodes[i].size;
    //     for(int j = 0; j < SIMD_WIDTH; ++j) {
    //         ispc_nodes[i].child[j] = nodes[i].child[j];
    //     }
    // }
    // ispc_ray.point = ispc_Vec3(ray.point);
    // ispc_ray.dir = ispc_Vec3(ray.dir);
    // ispc_ray.dist_bounds = ispc_Vec2(ray.dist_bounds);

            

    // for(int i = 0; i < SIMD_WIDTH; ++i) {
    //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
    //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
    // }

    // ispc::find_hit(0, ispc_ray, ispc_times, ispc_bbox, ispc_hits, ispc_triangles, ispc_ret, ispc_nodes);
    // Trace ret;
    // ret.hit = ispc_ret[0].hit;
    // ret.distance = ispc_ret[0].distance;
    // ret.position = Vec3(ispc_ret[0].position.x, ispc_ret[0].position.y, ispc_ret[0].position.z);
    // ret.normal = Vec3(ispc_ret[0].normal.x, ispc_ret[0].normal.y, ispc_ret[0].normal.z);
    // ret.origin = Vec3(ispc_ret[0].origin.x, ispc_ret[0].origin.y, ispc_ret[0].origin.z);
    // ret.material = ispc_ret[0].material;
    // return ret;

    auto is_leaf=[](Kmeans*kmeans){
        return kmeans->children[0]==NULL && kmeans->children[1]==NULL && kmeans->children[2]==NULL && kmeans->children[3]==NULL;
    };

    ispc::Ray ispc_ray;
    ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    // ispc::BBox *ispc_bbox = new ispc::BBox[nodes.size()];
    // // std::vector<ispc::BBox> ispc_bbox = new ispc::BBox[nodes.size()];
    // for(unsigned long i = 0; i < nodes.size(); i++) {
    //     ispc_bbox[i].min = ispc_Vec3(nodes[i].bbox.min);
    //     ispc_bbox[i].max = ispc_Vec3(nodes[i].bbox.max);
    // }
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2 *times = new Vec2[SIMD_WIDTH];

    // for(int i = 0; i < SIMD_WIDTH; ++i) {
    //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
    //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
    // }
            
    ispc_ray.point = ispc_Vec3(ray.point);
    ispc_ray.dir = ispc_Vec3(ray.dir);
    Trace ret;
    std::stack<Kmeans*> node_stack;
    node_stack.push(kmeans);
    while(!node_stack.empty()) {
        Kmeans* Kmeans_top = node_stack.top();
        node_stack.pop();
        // const Node4& node = nodes[node_idx];
        // printf("node idx: %zu\n", node_idx);
        // with early return
        // Vec2 times0{};
        // bool hit0 = node.bbox.hit(ray, times0);
        // if(!hit0 || (ret.hit && ret.distance <= times0.x))
        // {
        //     continue;
        // }
        // printf("here here111\n");
        if(is_leaf(Kmeans_top)) {
            // printf("is leaf is leaf here here111\n");
         
            for(size_t p=0;p<Kmeans_top->indexOfPrimitives.size();++p){
                Trace hit = primitives[Kmeans_top->indexOfPrimitives[p]].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            // printf("not leaf not leaf here here111\n");

            // Vec2 times[SIMD_WIDTH];
            // bool hits[SIMD_WIDTH];
            // for(int i = 0; i < SIMD_WIDTH; ++i) {
            //     hits[i] = nodes[node.child[i]].bbox.hit(ray, times[i]);
            // }
            // for(int i = 0; i < SIMD_WIDTH; ++i) {
            //     if(hits[i]) {
            //         if(node.child[i])
            //         {
            //             node_stack.push(node.child[i]);
            //         }
            //     }
            // }
            // printf("not leaf not leaf here here111\n");

            

            // for(int i = 0; i < SIMD_WIDTH; ++i) {
            //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
            //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
            // }

            // for(int i = 0; i < SIMD_WIDTH; ++i) {
            //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
            //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
            // }
            
            // std::vector<ispc::BBox> ispc_bbox4 = {ispc_bbox[node.child[0]], ispc_bbox[node.child[1]], ispc_bbox[node.child[2]], ispc_bbox[node.child[3]], ispc_bbox[node.child[4]], ispc_bbox[node.child[5]], ispc_bbox[node.child[6]], ispc_bbox[node.child[7]]};
            // printf("before ispc before ispc\n");

            // if(kmeans->children[j] != NULL && kmeans->cluster[j].indexOfPrimitives.size()) {
    //         visit = true;

    //         Vec2 times;
    //         bool hit = kmeans->cluster[j].representive.hit(ray, times);
    //         if(hit) {
    //             Trace hit = hit_kmeans(ray, kmeans->children[j]);
    //             ret = Trace::min(ret, hit);
    //         }
    //     }
            ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];
            std::vector<ispc::BBox> ispc_bbox4;
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(kmeans->children[i] && kmeans->cluster[i].indexOfPrimitives.size())
                {
                    ispc_bbox[i].min = ispc_Vec3(kmeans->cluster[i].representive.min);
                    ispc_bbox[i].max = ispc_Vec3(kmeans->cluster[i].representive.max);
                }
            }

            //  = {ispc_bbox[node.child[0]], ispc_bbox[node.child[1]], ispc_bbox[node.child[2]], ispc_bbox[node.child[3]]};
            ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);
            // printf("after ispc after ispc\n");
            
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(kmeans->children[i])
                    {
                        node_stack.push(kmeans->children[i]);
                    }
                }
            }




            // free(ispc_times);
            // free(ispc_bbox);
            // free(ispc_hits);
            // free(times);
            // printf("after free after free\n");
        }
    }
    return ret;

    //     for(size_t p=0;p<kmeans->cluster[i].indexOfPrimitives.size();++p){
    //       pTemp.push_back(kmeans->cluster[i].indexOfPrimitives[p]);
    //     kmeans->cluster[i].
    //     Trace hit = primitives[i].hit(ray);
    //     ret = Trace::min(ret, hit);

    //     kmeans->children=new Kmeans* [kmeans->m_K];
    // }
}

template<typename Primitive>
Trace BVH<Primitive>::find_hit(const Ray& ray, const size_t &idx) const {
    Trace ret;
    if (nodes[idx].is_leaf()) {
        for (unsigned long i = nodes[idx].start; i < nodes[idx].start + nodes[idx].size; i++) {
            Trace hit = primitives[i].hit(ray);
            ret = Trace::min(ret, hit);
        }
        return ret;
    } else {
        Vec2 times_l;
        Vec2 times_r;
        // pass reference
        bool hit_l = nodes[nodes[idx].l].bbox.hit(ray, times_l);
        bool hit_r = nodes[nodes[idx].r].bbox.hit(ray, times_r);

        if (hit_l && hit_r) {
            if (times_l.x < times_r.x) {
                Trace ret_l = find_hit(ray, nodes[idx].l);
                Trace ret_r = find_hit(ray, nodes[idx].r);
                return Trace::min(ret_l, ret_r);
            } else {
                Trace ret_r = find_hit(ray, nodes[idx].r);
                Trace ret_l = find_hit(ray, nodes[idx].l);
                return Trace::min(ret_l, ret_r);
            }
        } else if (hit_l) {
            return find_hit(ray, nodes[idx].l);
        } else if (hit_r) {
            return find_hit(ray, nodes[idx].r);
        } else {
            return ret;
        }
    }
}

template<typename Primitive>
Trace BVH<Primitive>::hit_queue(const Ray& ray) const {
    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node& node = nodes[node_idx];

        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x))
        {
            continue;
        }
        
        if(node.is_leaf()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                // printf("hit leaf node\n");
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            Vec2 times1, times2;
            bool hit1 = nodes[node.l].bbox.hit(ray, times1);
            bool hit2 = nodes[node.r].bbox.hit(ray, times2);
            if(hit1 && hit2) {
                size_t first = times1.x < times2.x ? node.l : node.r;
                size_t second = times1.x < times2.x ? node.r : node.l;
                node_stack.push(first);
                node_stack.push(second);
            } else if(hit1) {
                node_stack.push(node.l);
            } else if(hit2) {
                node_stack.push(node.r);
            }
        }
    }
    return ret;
}


template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive>
BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive>
bool BVH<Primitive>::Node::is_leaf() const {
    return l == r;
}


template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive>
BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive>
std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive>
void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

// template<typename Primitive>
// size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
//                                  const Mat4& trans) const {

//     std::stack<std::pair<size_t, size_t>> tstack;
//     tstack.push({root_idx, 0});
//     size_t max_level = 0;

//     if(nodes.empty()) return max_level;

//     while(!tstack.empty()) {

//         auto [idx, lvl] = tstack.top();
//         max_level = std::max(max_level, lvl);
//         const Node& node = nodes[idx];
//         tstack.pop();

//         Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
//         GL::Lines& add = lvl == level ? active : lines;

//         BBox box = node.bbox;
//         box.transform(trans);
//         Vec3 min = box.min, max = box.max;

//         auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

//         edge(min, Vec3{max.x, min.y, min.z});
//         edge(min, Vec3{min.x, max.y, min.z});
//         edge(min, Vec3{min.x, min.y, max.z});
//         edge(max, Vec3{min.x, max.y, max.z});
//         edge(max, Vec3{max.x, min.y, max.z});
//         edge(max, Vec3{max.x, max.y, min.z});
//         edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
//         edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
//         edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
//         edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
//         edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
//         edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

//         if(node.l && node.r) {
//             tstack.push({node.l, lvl + 1});
//             tstack.push({node.r, lvl + 1});
//         } else {
//             for(size_t i = node.start; i < node.start + node.size; i++) {
//                 size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
//                 max_level = std::max(c, max_level);
//             }
//         }
//     }
//     return max_level;
// }


// another visualize function for node4

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node4& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if((node.child[0] && (!node.child[1] && !node.child[2] && !node.child[3]))
        || (node.child[1] && (!node.child[0] && !node.child[2] && !node.child[3]))
        || (node.child[2] && (!node.child[0] && !node.child[1] && !node.child[3]))
        || (node.child[3] && (!node.child[0] && !node.child[1] && !node.child[2]))) {
            // std::cout << "Error: Node4 has only one child" << std::endl;
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c, max_level);
            }
        }
        else
        {
            for(size_t i = 0; i < SIMD_WIDTH; i++) {
                if(node.child[i]) {
                    tstack.push({node.child[i], lvl + 1});
                }
            }
        }
    }
    return max_level;
}


template<typename Primitive>
bool BVH<Primitive>::Node4::is_leaf4() const {
    return child[0] == child[1] && child[1] == child[2] && child[2] == child[3];
}

template<typename Primitive>
size_t BVH<Primitive>::new_node4(BBox box, size_t start, size_t size) {
    Node4 n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.child[0] = 0;
    n.child[1] = 0;
    n.child[2] = 0;
    n.child[3] = 0;
    nodes.push_back(n);
    return nodes.size() - 1;
}



template<typename Primitive>
Trace BVH<Primitive>::hit_queue4(const Ray& ray) const {
    // Trace ret;
    // std::stack<size_t> node_stack;
    // node_stack.push(0);
    // while(!node_stack.empty()) {
    //     size_t node_idx = node_stack.top();
    //     node_stack.pop();
    //     const Node4& node = nodes[node_idx];

    //     // with early return
    //     Vec2 times0{};
    //     bool hit0 = node.bbox.hit(ray, times0);
    //     if(!hit0 || (ret.hit && ret.distance <= times0.x))
    //     {
    //         continue;
    //     }
        
    //     if(node.is_leaf4()) {
    //         size_t node_end = node.start + node.size;
    //         for(size_t i = node.start; i < node_end; ++i) {
    //             // printf("hit leaf node\n");
    //             Trace hit = primitives[i].hit(ray);
    //             ret = Trace::min(ret, hit);
    //         }
    //     } else {
    //         Vec2 *times = new Vec2[4];
    //         bool *hit = new bool[4];
    //         for(int i = 0; i < 4; ++i) {
    //             hit[i] = nodes[node.child[i]].bbox.hit(ray, times[i]);
    //             if(hit[i]) {
    //                 node_stack.push(node.child[i]);
    //             }
    //         }
    //         // Vec2 times1, times2, times3, times4;
    //         // bool hit1 = nodes[node.l].bbox.hit(ray, times1);
    //         // bool hit2 = nodes[node.r].bbox.hit(ray, times2);
    //         // bool hit3 = nodes[node.l].bbox.hit(ray, times3);
    //         // bool hit4 = nodes[node.r].bbox.hit(ray, times4);
    //         // if(hit1 && hit2) {
    //         //     size_t first = times1.x < times2.x ? node.l : node.r;
    //         //     size_t second = times1.x < times2.x ? node.r : node.l;
    //         //     node_stack.push(first);
    //         //     node_stack.push(second);
    //         // } else if(hit1) {
    //         //     node_stack.push(node.l);
    //         // } else if(hit2) {
    //         //     node_stack.push(node.r);
    //         // }
    //     }
    // }
    // return ret;

    // auto ispc_Vec3 = [](const Vec3& v) {
	// 	ispc::Vec3 res;
	// 	res.x = v.x; res.y = v.y; res.z = v.z;
	// 	return res;
	// };

    // Trace ret;
    // std::stack<size_t> node_stack;
    // node_stack.push(0);
    // while(!node_stack.empty()) {
    //     size_t node_idx = node_stack.top();
    //     node_stack.pop();
    //     const Node4& node = nodes[node_idx];

    //     // with early return
    //     Vec2 times0{};
    //     bool hit0 = node.bbox.hit(ray, times0);
    //     if(!hit0 || (ret.hit && ret.distance <= times0.x))
    //     {
    //         continue;
    //     }
        
    //     if(node.is_leaf4()) {
    //         size_t node_end = node.start + node.size;
    //         for(size_t i = node.start; i < node_end; ++i) {
    //             // printf("hit leaf node\n");
    //             Trace hit = primitives[i].hit(ray);
    //             ret = Trace::min(ret, hit);
    //         }
    //     } else {
    //         ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    //         ispc::Ray ispc_ray;
    //         ispc_ray.point = ispc_Vec3(ray.point);
    //         ispc_ray.dir = ispc_Vec3(ray.dir);

    //         ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];

    //         for(int i = 0; i < SIMD_WIDTH; ++i) {
    //             ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
    //             ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
    //         }
            
    //         bool* ispc_hits = new bool[SIMD_WIDTH];
    //         ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);

    //         Vec2 *times = new Vec2[SIMD_WIDTH];
    //         for(int i = 0; i < SIMD_WIDTH; ++i) {
    //             times[i].x = ispc_times[i].x;
    //             times[i].y = ispc_times[i].y;
    //         }

    //         for(int i = 0; i < SIMD_WIDTH; ++i) {
    //             if(ispc_hits[i]) {
    //                 node_stack.push(node.child[i]);
    //             }
    //         }
    //     }
    // }
    // return ret;


    // // ispc::Ray ispc_ray;
    // // ispc_ray.point = ispc_Vec3(ray.point);
    // // ispc_ray.dir = ispc_Vec3(ray.dir);
    // // ispc::Node* ispc_nodes = new ispc::Node[nodes.size()];
    // // for(int i = 0; i < nodes.size(); i++) {
    // //     ispc_nodes[i].bbox.min = ispc_Vec3(nodes[i].bbox.min);
    // //     ispc_nodes[i].bbox.max = ispc_Vec3(nodes[i].bbox.max);
    // //     ispc_nodes[i].start = nodes[i].start;
    // //     ispc_nodes[i].size = nodes[i].size;
    // //     ispc_nodes[i].l = nodes[i].l;
    // //     ispc_nodes[i].r = nodes[i].r;
    // // }
    // // ispc::Trace ispc_ret;

    // // ispc::find_hit(ispc_ray, 0, ispc_nodes, ispc_ret);

    // // Trace ret;
    // // ret.hit = ispc_ret.hit;
    // // ret.distance = ispc_ret.distance;
    // // ret.position = Vec3(ispc_ret.position.x, ispc_ret.position.y, ispc_ret.position.z);
    // // ret.normal = Vec3(ispc_ret.normal.x, ispc_ret.normal.y, ispc_ret.normal.z);
    // // ret.origin = Vec3(ispc_ret.origin.x, ispc_ret.origin.y, ispc_ret.origin.z);
    // // ret.material = ispc_ret.material;
    // // return ret;






    auto ispc_Vec3 = [](const Vec3& v) {
		ispc::Vec3 res;
		res.x = v.x; res.y = v.y; res.z = v.z;
		return res;
	};

    // // for parallel whold tree
    // auto ispc_Vec2 = [](const Vec2& v) {
    //     ispc::Vec2 res;
    //     res.x = v.x; res.y = v.y;
    //     return res;
    // };

    // ispc::Ray ispc_ray;
    // ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    // ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];
    // bool* ispc_hits = new bool[SIMD_WIDTH];
    // Vec2 *times = new Vec2[SIMD_WIDTH];
    // ispc::Trace *ispc_ret = new ispc::Trace[1];
    // ispc::Triangle* ispc_triangles = new ispc::Triangle[primitives.size()];
    // ispc::Node* ispc_nodes = new ispc::Node[nodes.size()];
    // for(int i = 0; i < primitives.size(); i++) {
    //     ispc_triangles[i].v_0 = vertex_list[primitives[i].v0];
    //     ispc_triangles[i].v_1 = vertex_list[primitives[i].v1];
    //     ispc_triangles[i].v_2 = vertex_list[primitives[i].v2];
    // }
    // for(int i = 0; i < nodes.size(); i++) {
    //     ispc_nodes[i].bbox.min = ispc_Vec3(nodes[i].bbox.min);
    //     ispc_nodes[i].bbox.max = ispc_Vec3(nodes[i].bbox.max);
    //     ispc_nodes[i].start = nodes[i].start;
    //     ispc_nodes[i].size = nodes[i].size;
    //     for(int j = 0; j < SIMD_WIDTH; ++j) {
    //         ispc_nodes[i].child[j] = nodes[i].child[j];
    //     }
    // }
    // ispc_ray.point = ispc_Vec3(ray.point);
    // ispc_ray.dir = ispc_Vec3(ray.dir);
    // ispc_ray.dist_bounds = ispc_Vec2(ray.dist_bounds);

            

    // for(int i = 0; i < SIMD_WIDTH; ++i) {
    //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
    //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
    // }

    // ispc::find_hit(0, ispc_ray, ispc_times, ispc_bbox, ispc_hits, ispc_triangles, ispc_ret, ispc_nodes);
    // Trace ret;
    // ret.hit = ispc_ret[0].hit;
    // ret.distance = ispc_ret[0].distance;
    // ret.position = Vec3(ispc_ret[0].position.x, ispc_ret[0].position.y, ispc_ret[0].position.z);
    // ret.normal = Vec3(ispc_ret[0].normal.x, ispc_ret[0].normal.y, ispc_ret[0].normal.z);
    // ret.origin = Vec3(ispc_ret[0].origin.x, ispc_ret[0].origin.y, ispc_ret[0].origin.z);
    // ret.material = ispc_ret[0].material;
    // return ret;


    ispc::Ray ispc_ray;
    ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    ispc::BBox *ispc_bbox = new ispc::BBox[nodes.size()];
    for(unsigned long i = 0; i < nodes.size(); i++) {
        ispc_bbox[i].min = ispc_Vec3(nodes[i].bbox.min);
        ispc_bbox[i].max = ispc_Vec3(nodes[i].bbox.max);
    }
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2 *times = new Vec2[SIMD_WIDTH];

            
    ispc_ray.point = ispc_Vec3(ray.point);
    ispc_ray.dir = ispc_Vec3(ray.dir);
    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node4& node = nodes[node_idx];
        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x))
        {
            continue;
        }
        if(node.is_leaf4()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            std::vector<ispc::BBox> ispc_bbox4;
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(node.child[i])
                {
                    ispc_bbox4.push_back(ispc_bbox[node.child[i]]);
                }
            }
            ispc::bbox_hit(ispc_ray, ispc_bbox4.data(), ispc_times, ispc_hits);
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(node.child[i])
                    {
                        node_stack.push(node.child[i]);
                    }
                }
            }

        }
    }
    return ret;
}


template<typename Primitive>
void BVH<Primitive>::SAH4(const size_t &idx, const size_t &max_leaf_size) {
    if (nodes[idx].size <= max_leaf_size) {
        return;
    }

    // Create bounding boxes for children
    BBox *split_bbox = new BBox[SIMD_WIDTH];
    size_t *rangel = new size_t[SIMD_WIDTH];
    size_t *ranger = new size_t[SIMD_WIDTH];
    size_t *node_addr = new size_t[SIMD_WIDTH];
    for(int i = 0; i < SIMD_WIDTH; ++i) {
        node_addr[i] = new_node4();
    }

    size_t *startl = new size_t[SIMD_WIDTH];
    size_t *startr = new size_t[SIMD_WIDTH];
    bucket_split(idx, rangel[0], ranger[0], split_bbox[0], split_bbox[1]);

    startl[0] = nodes[idx].start;
    startr[0] = startl[0] + rangel[0];

    for(int i = 0; i < 2; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i%2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i/2))];
            nodes[node_addr[i]].size = rangel[int(floor(i/2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i/2))];
            nodes[node_addr[i]].size = ranger[int(floor(i/2))];
        }
    }


    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    for(int i = 0; i < 2; ++i) {
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }
    for(int i = 0; i < 4; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i%2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i/2))];
            nodes[node_addr[i]].size = rangel[int(floor(i/2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i/2))];
            nodes[node_addr[i]].size = ranger[int(floor(i/2))];
        }
    }

    for(int i = 0; i < 4; ++i) {
        if(nodes[node_addr[i]].size)
        {
            nodes[idx].child[i] = node_addr[i];
        }
    }

    for(int i = 0; i < 4; ++i) {
        if(nodes[node_addr[i]].size)
        {
            SAH4(nodes[idx].child[i], max_leaf_size);
        }
    }
}






template<typename Primitive>
bool BVH<Primitive>::Node8::is_leaf8() const {
    return child[0] == 0 && child[0] == child[1] && child[1] == child[2] && child[2] == child[3] && child[3] == child[4] && child[4] == child[5] && child[5] == child[6] && child[6] == child[7];
}

template<typename Primitive>
size_t BVH<Primitive>::new_node8(BBox box, size_t start, size_t size) {
    Node8 n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    for(size_t i = 0; i < SIMD_WIDTH; i++) {
        n.child[i] = 0;
    }
    nodes.push_back(n);
    return nodes.size() - 1;
}



template<typename Primitive>
Trace BVH<Primitive>::hit_queue8(const Ray& ray) const {
    auto ispc_Vec3 = [](const Vec3& v) {
		ispc::Vec3 res;
		res.x = v.x; res.y = v.y; res.z = v.z;
		return res;
	};

    ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2 *times = new Vec2[SIMD_WIDTH];
    ispc::Ray ispc_ray;

    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node8& node = nodes[node_idx];
        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x))
        {
            continue;
        }
        if(node.is_leaf8()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            ispc_ray.point = ispc_Vec3(ray.point);
            ispc_ray.dir = ispc_Vec3(ray.dir);

            

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
                ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
            }
            
            ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(node.child[i])
                    {
                        node_stack.push(node.child[i]);
                    }
                }
            }
        }
    }
    return ret;
}


template<typename Primitive>
void BVH<Primitive>::SAH8(const size_t & idx, const size_t & max_leaf_size) {
    if (nodes[idx].size <= max_leaf_size) {
        return;
    }
    
    // Create bounding boxes for children
    BBox *split_bbox = new BBox[SIMD_WIDTH];
    size_t *rangel = new size_t[SIMD_WIDTH];
    size_t *ranger = new size_t[SIMD_WIDTH];
    size_t *node_addr = new size_t[SIMD_WIDTH];
    for(int i = 0; i < SIMD_WIDTH; ++i) {
        node_addr[i] = new_node8();
    }

    size_t *startl = new size_t[SIMD_WIDTH];
    size_t *startr = new size_t[SIMD_WIDTH];
    bucket_split(idx, rangel[0], ranger[0], split_bbox[0], split_bbox[1]);

    startl[0] = nodes[idx].start;
    startr[0] = startl[0] + rangel[0];

    for(int i = 0; i < 2; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i%2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i/2))];
            nodes[node_addr[i]].size = rangel[int(floor(i/2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i/2))];
            nodes[node_addr[i]].size = ranger[int(floor(i/2))];
        }
    }


    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    for(int i = 0; i < 2; ++i) {
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }
    for(int i = 0; i < 4; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i%2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i/2))];
            nodes[node_addr[i]].size = rangel[int(floor(i/2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i/2))];
            nodes[node_addr[i]].size = ranger[int(floor(i/2))];
        }
    }


    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    bucket_split(node_addr[2], rangel[2], ranger[2], split_bbox[4], split_bbox[5]);
    bucket_split(node_addr[3], rangel[3], ranger[3], split_bbox[6], split_bbox[7]);
    for(int i = 0; i < 4; ++i) {
        // printf("here %d\n", i);
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }

    // printf("reach here\n");
    for(int i = 0; i < 8; ++i) {
        // printf("node_addr[i] = %zu\n", node_addr[i]);
        if (node_addr[i]<0xffffffff)
        {
            nodes[node_addr[i]].bbox = split_bbox[i];
            if(i%2 == 0) {
                nodes[node_addr[i]].start = startl[int(floor(i/2))];
                nodes[node_addr[i]].size = rangel[int(floor(i/2))];
            } else {
                nodes[node_addr[i]].start = startr[int(floor(i/2))];
                nodes[node_addr[i]].size = ranger[int(floor(i/2))];
            }
        }
    }

    for(int i = 0; i < 8; ++i) {
        if(nodes[node_addr[i]].size)
        {
            nodes[idx].child[i] = node_addr[i];
        }
    }

    for(int i = 0; i < 8; ++i) {
        if(nodes[node_addr[i]].size)
        {
            SAH8(nodes[idx].child[i], max_leaf_size);
        }
    }
}


template<typename Primitive>
void BVH<Primitive>::bucket_split(const size_t& idx, size_t& rangel, size_t& ranger, BBox& split_leftBox, BBox& split_rightBox) {
    int bucket_num = 8;
    BBox bbox = nodes[idx].bbox;
    int min_cost_axis = 0;
    float min_cost = 0x7fffffff;
    float min_cost_split = 0;
    // can't use buckets[3][bucket_num], will cause bus fault
    for (int i = 0; i < 3; i++) {
        BBox buckets[bucket_num];
        int prim_count[bucket_num];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < bucket_num; j++) {
                prim_count[j] = 0;
            }
        }
        float my_min = bbox.min[i];
        float my_max = bbox.max[i];
        float interval = (my_max - my_min) / (bucket_num + 0.0);
        for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
            Primitive& p = primitives[j];
            BBox pbb = p.bbox();
            
            int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
            bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);
            buckets[bucket_idx].enclose(pbb);
            prim_count[bucket_idx]++;
        }
        
        
        for (int j = 0; j < bucket_num - 1; j++) {
            float left_surface_area = 0;
            float right_surface_area = 0;
            float left_prim_count = 0;
            float right_prim_count = 0;

            BBox left_bbox;
            BBox right_bbox;
    
            for (int k = 0; k <= j; k++) {
                left_bbox.enclose(buckets[k]);
                left_prim_count += prim_count[k];
            }
            for (int k = j + 1; k < bucket_num; k++) {
                right_bbox.enclose(buckets[k]);
                right_prim_count += prim_count[k];
            }

            left_surface_area = left_bbox.surface_area();
            right_surface_area = right_bbox.surface_area();

            float total_cost = left_surface_area * left_prim_count + right_surface_area * right_prim_count;
            if (total_cost < min_cost) {
                min_cost = total_cost;
                min_cost_axis = i;
                min_cost_split = my_min + (j + 1) * interval;
                split_leftBox = left_bbox;
                split_rightBox = right_bbox;
                rangel = left_prim_count;
                ranger = right_prim_count;
            }
        }
    }

    // need to reorganize primitives so that the children are contiguous ranges of primitives
    int first = nodes[idx].start;
    for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) 
    {
        Primitive& p = primitives[j];
        BBox pbb = p.bbox();
        if (pbb.center()[min_cost_axis] < min_cost_split)
        {
            std::swap(primitives[j], primitives[first]);
            ++first;
        }
    }

}

} // namespace PT
