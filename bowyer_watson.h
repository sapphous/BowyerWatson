
#ifndef BOWYER_WATSON_H_
#define BOWYER_WATSON_H_


#include<vector>
#include<map>
#include<Eigen/Dense>

//BowyerWatson control class WARNING: make all methods static, or perhaps replace with a namespace? No real point to having a class here except for convenience.
struct BowyerWatson
{
    class Triangle;
    struct PointComparator;
    struct Mesh;

    static const unsigned int cache_line_size = 64;
    static const unsigned int point_size = sizeof(std::pair<Eigen::Vector3f, unsigned int>);

    //generate how_many points on a sphere, placing them into empty_vector (the unsigned ints are for storing the triangle to which they belong later)
    void generate_points(unsigned int how_many, std::vector<std::pair<Eigen::Vector3f, unsigned int> > &empty_vector) const;

    //create an initial triangulation (8 triangles comprising two square pyramids sharing a base, using points from the pointset as vertices)
    //returns a pair consisting of the first and last triangle in a linked list of 8 triangles; the points passed are sorted and contained within them
    std::pair<Triangle*, Triangle*> create_initial_triangles(std::vector<std::pair<Eigen::Vector3f, unsigned int> > &points) const;
    std::pair<Triangle*, Triangle*> create_initial_triangles(unsigned int how_many) const; //calls generate_points internally

    //create an initial triangulation (8 triangles comprising two square pyramids sharing a base, using the intersections of the x,y,z axes with the unit sphere as the vertices)
    //returns a pair consisting of the first and last triangle in a linked list of 8 triangles; the points passed are sorted and contained within them
    std::pair<Triangle*, Triangle*> create_cardinal_triangles(std::vector<std::pair<Eigen::Vector3f, unsigned int> > &points) const;
    std::pair<Triangle*, Triangle*> create_cardinal_triangles(unsigned int how_many_minus_six) const; //calls generate_points internally

    //perform the Bowyer-Watson algorithm, returning a pointer to the first triangle in the resulting null-terminated linked list of triangles
    Triangle *perform(Triangle *start_triangle, Triangle *last_triangle);
    Triangle *perform(std::pair<Triangle*, Triangle*> pair)                                 {   return(perform(pair.first, pair.second));   }

    //get a Mesh from the triangles, optionally deleting them afterwards WARNING: total_points is currently unused, but can be used for optimization
    Mesh get_mesh(Triangle *first, bool destructive = false, unsigned int total_points = 0);

    //takes the output of perform() and constructs a new pointset by replacing each point with the weighted average of the centroids of all triangles of which it is a vertex
    void relax_points(Triangle *start_triangle, std::vector<std::pair<Eigen::Vector3f, unsigned int> > &empty_vector, bool kill_triangles = true) const;
};



//class representing a spherical triangle on the unit sphere
class BowyerWatson::Triangle
{
    std::pair<Eigen::Vector3f, float> calculate_circumsphere() const;
public:
    Eigen::Vector3f vertices[3];    //vertices in counterclockwise order
    Triangle *neighbors[3];         //neighbors[n] is along the edge given by vertices[(n+1)%3] and vertices[(n+2)%3]

    std::pair<Eigen::Vector3f, float> circumsphere; //(circumcenter, radius^2 of circumsphere)

    unsigned int pointcount;        //the number of unprocessed points inside this triangle
    std::pair<Eigen::Vector3f, unsigned int> *points;       //an array of unprocessed points inside this triangle; the unsigned int will store the index of the new triangle for the point to be moved into as the algorithm progresses
    Triangle *prev;                 //the previous triangle (all triangles are kept track of in a linked list)
    Triangle *next;                 //the next triangle (all triangles are kept track of in a linked list)

    bool tagged;                    //whether we have been tagged as processed during a neighbor check
    bool is_bad;                    //whether we are a bad triangle or not in the current neighbor check (bad triangles will be removed after sorting their points to new ones)

    //arguments: vertices 0, 1, and 2 in counter-clockwise order (viewed from above, outside the unit sphere); n0 is the neighbor triangle along edge (v1,v2), which is convenient within the algorithm to initialize in the constructor: can be initialized to nullptr otherwise
    Triangle(const Eigen::Vector3f &v0, const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, Triangle *n0);
    ~Triangle(); //frees points array but does not screw with prev and next or do anything else

    //begins a neighbor_check, calling neighbor_check on neighboring triangles recursively; a neighbor-check finds all triangles whose circumspheres include origin;
    //  total_points should start at 0, and will be incremented to count the points;
    //  bads should begin empty and will be filled with triangles whose circumspheres include origin (these will be destroyed);
    //  polygon should begin empty and will be filled with (v, (t,i)) pairs representing the non-bad triangles bordering the star-shaped hole resulting when the bad triangles are removed:
    //    t is a pointer to the non-bad triangle
    //    i is the index of the bad neighbor (vertices (i+1)%3 and (i+2)%3 define the edge shared with the bad neighbor)
    //    v is the point on the edge shared with the bad neighbor which is immediately adjacent in the counter-clockwise direction to the other point on that edge, from the perspective of the hole; thus it is immediately adjacent in the clockwise direction from the perspective of t, i.e. it is vertex (i+2)%3
    void start_neighbor_check(Eigen::Vector3f &origin, unsigned int &total_points, std::vector<Triangle*> &bads, std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator> &polygon);

    //recursively called by itself and start_neighbor_check
    void neighbor_check(Eigen::Vector3f &point, unsigned int &total_points, std::vector<Triangle*> &bads, std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator> &polygon);

    //returns true if and only if point is contained by this triangle; used by create_initial_triangles() and useful as a general utility, but not used by the BW algorithm itself
    bool contains(const Eigen::Vector3f &point) const;
};


//a comparison class to allow Eigen::Vector3f objects to be used in STL containers
struct BowyerWatson::PointComparator
{
    bool operator()(const Eigen::Vector3f &a, const Eigen::Vector3f &b)
    {
        if(a[0] != b[0])
            return(a[0] < b[0]);
        if(a[1] != b[1])
            return(a[1] < b[1]);
        return(a[2] < b[2]);
    }
};

//a triangular mesh
struct BowyerWatson::Mesh
{
    std::vector<float> vertices;            //(vertices[3n], vertices[3n+1], vertices[3n+2]) = the nth point's x,y,z coordinates
    std::vector<unsigned int> triangles;    //( triangles[3n]th point, triangles[3n+1]th point, triangles[3n+2]th point ) = the nth triangle's vertices
};




#endif /* BOWYER_WATSON_H_ */
