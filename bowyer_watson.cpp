
#include<stdlib.h>
#include "bowyer_watson.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

//struct for Intel TBB to iterate over points in bad triangles (i.e. triangles to be removed from the triangulation)
struct PointInterval
{
    //MOOSE WARNING: this 7 is a result of some empirical tests & is tweakable. also points_per_run*sizeof(...) needs to divide cache_line_size to avoid cache line conflicts
    static const unsigned int points_per_run = 7 * BowyerWatson::cache_line_size / sizeof(std::pair<Eigen::Vector3f, unsigned int>);
    static bool worth_parallelizing(unsigned int total_points);

    struct CountNode //binary tree of unsigned-int-bearing nodes constructed by intel TBB before thead fork; each thread can count points in its own variable, and after they're done the main thread will sum them up and delete the structure
    {
        unsigned int count;
        CountNode *children[2];
        CountNode() : count(0), children{nullptr, nullptr}  {}
        ~CountNode()                    {   if(children[0] != nullptr) delete children[0]; if(children[1] != nullptr) delete children[1];   }
        void sum_into(unsigned int &v)  {   v += count; if(children[0] != nullptr) children[0]->sum_into(v); if(children[1] != nullptr) children[1]->sum_into(v);   }
    };

    CountNode *countnode;
    std::vector<BowyerWatson::Triangle*> &sources; //vector of triangles from which to take points
    unsigned int total_points;      //how many points are in this interval?
    unsigned int triangle_start;    //at what triangle in the sources vector does this interval start?
    unsigned int triangle_end;      //at what triangle in the sources vector does this interval end? (inclusive)
    unsigned int point_start;       //at what index in the point array within triangle_start does the interval start?
    unsigned int point_end;         //at what index in the point array within triangle_end does the interval end? (exclusive)



    bool empty() const                              {   return(total_points == 0);   }
    bool is_divisible() const                       {   return(total_points > 2 * points_per_run);   }

    PointInterval(PointInterval& r, tbb::split) :
        countnode(new CountNode),
        sources(r.sources), total_points(r.total_points),
        triangle_start(r.triangle_start), triangle_end(r.triangle_end),
        point_start(r.point_start), point_end(r.point_end)
    {
        r.countnode->children[1] = countnode;
        r.countnode = (r.countnode->children[0] = new CountNode);
        while(total_points > r.total_points / 2)
        {
            if(sources[triangle_start]->pointcount - point_start <= points_per_run)
            {
                total_points -= sources[triangle_start]->pointcount - point_start;
                ++triangle_start; //we never need to worry about this taking us past sources.back(), since we're only going halfway through total_points
                point_start = 0;
            }
            else
            {
                total_points -= points_per_run;
                point_start += points_per_run;
            }
        }
        r.total_points -= total_points;
        r.triangle_end = triangle_start;
        r.point_end = point_start;
    }
    PointInterval(std::vector<BowyerWatson::Triangle*> &sources_, CountNode *countnode_) : //sources_ is (bad) triangles to copy points from, countnode_ is Countnode to count points in
        countnode(countnode_),
        sources(sources_), total_points(0),
        triangle_start(0), triangle_end(sources.size() - 1),
        point_start(0), point_end(sources.back()->pointcount)
    {
        for(unsigned int i = 0; i < sources.size(); ++i)
            total_points += sources[i]->pointcount;
    }
};

bool PointInterval::worth_parallelizing(unsigned int total_points)
{
    return(total_points > 4 * points_per_run); //MOOSE WARNING: this is just a rough empirically good choice
}


//struct for finding those points within a PointInterval contained within a given triangle
//WARNING: only to be used within BW algorithm, i.e. in perform(); makes the assumption that we are dealing with a triangular slice of a star-shaped region in determining containment
struct PointSorter
{
    unsigned int tag;       //the tag to put in the points array to mark a point as belonging to our triangle
    Eigen::Vector3f n[2];   //the normal vectors to the appropriate edges of our triangle (we know we are in a star-shaped region, so the edge opposite to our center need not be checked)

    PointSorter(unsigned int tag_, BowyerWatson::Triangle *container) :
        tag(tag_), n{
            container->vertices[0].cross(container->vertices[1]),
            container->vertices[2].cross(container->vertices[0])
        }
    {}
    //takes the points in a PointInteval and checks whether they are within the triangle passed to the constructor
    void operator()(const PointInterval &interval) const
    {
        unsigned int count(0);
        std::pair<Eigen::Vector3f, unsigned int> *arr;
        unsigned int pindex = interval.point_start;
        unsigned int stoppoint;
        for(unsigned int trindex = interval.triangle_start; trindex <= interval.triangle_end; ++trindex)
        {
            //determine how many points to iterate
            if(trindex == interval.triangle_end)
                stoppoint = interval.point_end;
            else
                stoppoint = interval.sources[trindex]->pointcount;
            //create a pointer to the next point to read, and begin reading
            arr = interval.sources[trindex]->points;
            for(; pindex < stoppoint; ++pindex)
            {
                if(arr[pindex].second == (unsigned int)(-1) && arr[pindex].first.dot(n[0]) > 0 && arr[pindex].first.dot(n[1]) > 0)
                {
                    arr[pindex].second = tag;
                    ++count;
                }
            }
            //reset our point index to 0
            pindex = 0;
        }
        interval.countnode->count = count;
    }
};

//struct for moving points from a set of to-be-destroyed triangles into the appropriate new triangles
struct PointMover
{
    std::vector<BowyerWatson::Triangle*> &sources;
    std::vector<BowyerWatson::Triangle*> &destinations;

    PointMover(std::vector<BowyerWatson::Triangle*> &sources_, std::vector<BowyerWatson::Triangle*> &destinations_) :
        sources(sources_), destinations(destinations_)
    {}

    void operator()(const tbb::blocked_range<unsigned int> &interval) const
    {
        //create fun local variables for possible efficiency
        unsigned int begin = interval.begin();
        unsigned int end = interval.end();
        //allocate destination arrays and create a vector of pointers to their first elements (these will be incremented as we insert)
        std::vector<std::pair<Eigen::Vector3f, unsigned int>*> dest_ptr; //dest_ptr[x-begin] is the address where the next point inside destinations[x] should be copied
        for(unsigned int k = begin; k < end; ++k)
        {
            //align array to cache boundary; the 2nd argument is size, & to ensure it is divisible by alignment (for some godawful reason the standard requires this),
            //we perform Ceiling(needed bytes / cache line size) * (cache line size)
            dest_ptr.push_back(destinations[k]->points = (std::pair<Eigen::Vector3f, unsigned int>*) aligned_alloc(
                BowyerWatson::cache_line_size,
                ((destinations[k]->pointcount * BowyerWatson::point_size + BowyerWatson::cache_line_size - 1) / BowyerWatson::cache_line_size) * BowyerWatson::cache_line_size
            ));
        }
        //copy points in sources[n] for all n whenever the point belongs in one of our destinations (those between destinations[begin] and destinations[end])
        std::pair<Eigen::Vector3f, unsigned int> *to_check;
        for(unsigned int i = 0; i < sources.size(); ++i)
        {
            to_check = sources[i]->points;
            for(unsigned int j = 0; j < sources[i]->pointcount; ++to_check, ++j)
                if(to_check->second >= begin && to_check->second < end)
                    *(dest_ptr[to_check->second - begin]++) = std::make_pair(to_check->first, (unsigned int)(-1));
        }
    }
};


std::pair<Eigen::Vector3f, float> BowyerWatson::Triangle::calculate_circumsphere() const
{
    Eigen::Vector3f ac = vertices[2] - vertices[0];
    Eigen::Vector3f ab = vertices[1] - vertices[0];
    Eigen::Vector3f abXac = ab.cross(ac);
    Eigen::Vector3f aToCenter = (abXac.cross(ab) * ac.dot(ac) + ac.cross(abXac) * ab.dot(ab)) / (2.0f * abXac.dot(abXac));
    return(std::make_pair(
        vertices[0] + aToCenter,
        aToCenter.dot(aToCenter)
    ));
}


BowyerWatson::Triangle::Triangle(const Eigen::Vector3f &v0, const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, Triangle *n0) :
    vertices{v0, v1, v2}, neighbors{n0, nullptr, nullptr},
    circumsphere(calculate_circumsphere()),
    pointcount(0), points(nullptr), prev(nullptr), next(nullptr),
    tagged(false), is_bad(false)
{
}

BowyerWatson::Triangle::~Triangle()
{
    if(points != nullptr)
        free(points);
}


void BowyerWatson::Triangle::start_neighbor_check(Eigen::Vector3f &origin, unsigned int &total_points, std::vector<Triangle*> &bads, std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator> &polygon)
{
        //we're bad
        tagged = true;
        is_bad = true;
        bads.push_back(this);
        total_points += pointcount;

        neighbors[0]->neighbor_check(origin, total_points, bads, polygon);
        if(!neighbors[1]->tagged)
            neighbors[1]->neighbor_check(origin, total_points, bads, polygon);
        if(!neighbors[2]->tagged)
            neighbors[2]->neighbor_check(origin, total_points, bads, polygon);
}

void BowyerWatson::Triangle::neighbor_check(Eigen::Vector3f &origin, unsigned int &total_points, std::vector<Triangle*> &bads, std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator> &polygon)
{
    //check whether we're bad
    Eigen::Vector3f diff = origin - circumsphere.first;
    if(diff.dot(diff) >= circumsphere.second)
    {
        //we're not bad
        tagged = true;
        is_bad = false;
        //if a neighbor is already known to be bad, it is our responsibility to put our border with that neighbor into polygon
        if(neighbors[0]->tagged && neighbors[0]->is_bad)            polygon.insert(std::make_pair(vertices[2], std::make_pair(this, 0)));
        if(neighbors[1]->tagged && neighbors[1]->is_bad)            polygon.insert(std::make_pair(vertices[0], std::make_pair(this, 1)));
        if(neighbors[2]->tagged && neighbors[2]->is_bad)            polygon.insert(std::make_pair(vertices[1], std::make_pair(this, 2)));
        //no recursive calls to neighbor_check; all bad triangles have a bad neighbor, so there is no need to recurse from a non-bad triangle
    }
    else
    {
        //we're bad
        tagged = true;
        is_bad = true;
        bads.push_back(this);
        total_points += pointcount;
        //if a neighbor is already known to be non-bad, it is our responsibility to put our border with that neighbor into polygon
        if(neighbors[0]->tagged && !neighbors[0]->is_bad)           polygon.insert(std::make_pair(vertices[1], std::make_pair(neighbors[0],
                                                                        neighbors[0]->vertices[0] == vertices[1] ? 1 :
                                                                        neighbors[0]->vertices[1] == vertices[1] ? 2 : 0
                                                                    )));
        if(neighbors[1]->tagged && !neighbors[1]->is_bad)           polygon.insert(std::make_pair(vertices[2], std::make_pair(neighbors[1],
                                                                        neighbors[1]->vertices[0] == vertices[2] ? 1 :
                                                                        neighbors[1]->vertices[1] == vertices[2] ? 2 : 0
                                                                    )));
        if(neighbors[2]->tagged && !neighbors[2]->is_bad)           polygon.insert(std::make_pair(vertices[0], std::make_pair(neighbors[2],
                                                                        neighbors[2]->vertices[0] == vertices[0] ? 1 :
                                                                        neighbors[2]->vertices[1] == vertices[0] ? 2 : 0
                                                                    )));
        //recursively call neighbor_check on neighbors for whom it hasn't yet been called
        //WARNING: must check tagged before each of these calls, because the call to neighbors[0]->neighbor_check might tag neighbors[1], etc.
        if(!neighbors[0]->tagged)                                   neighbors[0]->neighbor_check(origin, total_points, bads, polygon);
        if(!neighbors[1]->tagged)                                   neighbors[1]->neighbor_check(origin, total_points, bads, polygon);
        if(!neighbors[2]->tagged)                                   neighbors[2]->neighbor_check(origin, total_points, bads, polygon);
    }
}

bool BowyerWatson::Triangle::contains(const Eigen::Vector3f &point) const
{
    Eigen::Matrix3f m;
    //check plane including vertices 0 and 1
    m.row(0) = vertices[0];
    m.row(1) = vertices[1];
    m.row(2) = point;
    if(m.determinant() < 0)
        return(false);
    //check plane including vertices 1 and 2
    m.row(0) = vertices[1];
    m.row(1) = vertices[2];
    m.row(2) = point;
    if(m.determinant() < 0)
        return(false);
    //check plane including vertices 2 and 0
    m.row(0) = vertices[2];
    m.row(1) = vertices[0];
    m.row(2) = point;
    return(m.determinant() >= 0);
}


BowyerWatson::Triangle *BowyerWatson::perform(Triangle *start_triangle, Triangle *last_triangle)
{
    //initial setup
    Triangle *triangle = start_triangle;
    Eigen::Vector3f point;
    std::vector<Triangle*> bads;
    std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator> polygon;    //star-shaped border segments parameterized by clockwise-most vertex of each
    std::vector<Triangle*> created;
    unsigned int total_points;
    //iterate through all triangles in the linked list
    while(triangle != 0)
    {
        //skip triangle if it contains no points
        if(triangle->pointcount == 0)
        {
            triangle = triangle->next;
            continue;
        }
        point = triangle->points[triangle->pointcount - 1].first;   //get the last point in the array to be our new vertex
        --triangle->pointcount;                                     //array is deleted using free, so this variable can lie; make it lie to hide the point we're using as a vertex so that it won't be copied into newly created triangles
        //find bad triangles and bounding polygon
        total_points = 0;
        bads.clear();
        polygon.clear();
        triangle->start_neighbor_check(point, total_points, bads, polygon);
        //create new triangles
        created.clear();
        std::map<Eigen::Vector3f, std::pair<Triangle*, unsigned char>, PointComparator>::iterator cur_edge = polygon.begin();
        do
        {
            //untag the triangle that won't get deleted
            cur_edge->second.first->tagged = false;
            //create the new triangle
            created.push_back(cur_edge->second.first->neighbors[cur_edge->second.second] = new Triangle(
                point,
                cur_edge->second.first->vertices[cur_edge->second.second > 0 ? cur_edge->second.second - 1 : cur_edge->second.second + 2],
                cur_edge->second.first->vertices[cur_edge->second.second > 1 ? cur_edge->second.second - 2 : cur_edge->second.second + 1],
                cur_edge->second.first
            ));
            //move counterclockwise around the perimeter of our polygon
            cur_edge = polygon.find(cur_edge->second.first->vertices[cur_edge->second.second > 1 ? cur_edge->second.second - 2 : cur_edge->second.second + 1]);
        } while(cur_edge != polygon.begin()); //once we've travelled all the way around the polygon back to our starting point, we stop
        //connect new triangles to one another
        created[0]->neighbors[1] = created[1];
        created[0]->neighbors[2] = created.back();
        created.back()->neighbors[1] = created[0];
        created.back()->neighbors[2] = created[created.size()-2]; //there are always at least 3 triangles created, since the point lies in *triangle
        for(unsigned int i = 1; i < created.size() - 1; ++i)
        {
            created[i]->neighbors[1] = created[i+1];
            created[i]->neighbors[2] = created[i-1];
        }
        //determine which new triangles to put points into
        if(PointInterval::worth_parallelizing(total_points))
        {
            for(unsigned int i = 0; i < created.size(); ++i)
            {
                //go through all points in the bad triangles to find those that should be in created[i]
                PointInterval::CountNode *countnode = new PointInterval::CountNode();
                tbb::parallel_for(PointInterval(bads, countnode), PointSorter(i, created[i]));
                //record the count of those points
                countnode->sum_into(created[i]->pointcount);
                delete countnode;
            }
        }
        else
        {
            for(unsigned int i = 0; i < created.size(); ++i)
            {
                //go through all points in the bad triangles to find those that should be in created[i]
                PointInterval::CountNode *countnode = new PointInterval::CountNode();
                PointSorter(i, created[i])(PointInterval(bads, countnode));
                //record the count of those points
                countnode->sum_into(created[i]->pointcount);
                delete countnode;
            }
        }
        //copy points into new triangles MOOSE WARNING: call a function to check whether parallelization is worth it, rather than doing it inline (450 is empirically tested as a pretty good cutoff; may not be optimal)
        if(total_points > 450)
            tbb::parallel_for(tbb::blocked_range<unsigned int>(0, created.size()), PointMover(bads, created));
        else
            PointMover(bads, created)(tbb::blocked_range<unsigned int>(0, created.size()));
        //add new triangles to the linked list (empty ones to the start, others to the end)
        for(unsigned int i = 0; i < created.size(); ++i)
        {
            if(created[i]->pointcount == 0)
            {
                start_triangle->prev = created[i];
                created[i]->next = start_triangle;
                start_triangle = created[i];
            }
            else
            {
                last_triangle->next = created[i];
                created[i]->prev = last_triangle;
                last_triangle = created[i];
            }
        }
        //destroy old bad triangles
        for(unsigned int i = 0; i < bads.size(); ++i)
        {
            if(bads[i] == start_triangle)   start_triangle = bads[i]->next;
            else                            bads[i]->prev->next = bads[i]->next;    //prev != nullptr since we are not *start_triangle
            if(bads[i] == last_triangle)    last_triangle = bads[i]->prev;
            else                            bads[i]->next->prev = bads[i]->prev;    //next != nullptr since we are not *last_triangle
            if(bads[i] == triangle)
                triangle = triangle->next;
            delete bads[i];
        }
    }
    return(start_triangle);
}


void BowyerWatson::generate_points(unsigned int how_many, std::vector<std::pair<Eigen::Vector3f, unsigned int> > &created) const
{
    created.clear();
    created.reserve(how_many);
    float rand1;
    float rand2;
    for(unsigned int i = 0; i < how_many; ++i)
    {
        //create the new point
        rand1 = ((float)rand())/((float)RAND_MAX) * 2.0f - 1.0f;
        rand2 = ((float)rand())/((float)RAND_MAX) * 3.14159265358979323846f * 2.0f;
        created.push_back(std::make_pair(Eigen::Vector3f(
            sqrtf(1.0f - rand1*rand1) * cos(rand2),
            sqrtf(1.0f - rand1*rand1) * sin(rand2),
            rand1
        ), 0));
        //guess which triangle contains the new point under the assumption that the critical points end up exactly at (1,0,0), (-1,0,0), (0,1,0), etc.
        //WARNING: the indices here refer to the by_number array of triangles, created later
        if(created.back().first[0] < 0) created.back().second |= 1;
        if(created.back().first[1] < 0) created.back().second |= 2;
        if(created.back().first[2] < 0) created.back().second |= 4;
    }
}

std::pair<BowyerWatson::Triangle*, BowyerWatson::Triangle*> BowyerWatson::create_initial_triangles(unsigned int how_many) const
{
    std::vector<std::pair<Eigen::Vector3f, unsigned int> > created;
    generate_points(how_many, created);
    return(create_initial_triangles(created));
}

std::pair<BowyerWatson::Triangle*, BowyerWatson::Triangle*> BowyerWatson::create_initial_triangles(std::vector<std::pair<Eigen::Vector3f, unsigned int> > &created) const
{
    unsigned int how_many = created.size();
    //find the critical points closest to the intersections of the x,y,z axes with the unit sphere
    float max[3] = { -2.0, -2.0, -2.0 };
    float min[3] = { 2.0, 2.0, 2.0 };
    unsigned int critical_index[6] = { 0, 0, 0, 0, 0, 0 };
    unsigned int *max_index = &critical_index[0];
    unsigned int *min_index = &critical_index[3];
    for(unsigned int i = 0; i < how_many; ++i)
    {
        //update the critical points
        for(unsigned int d = 0; d < 3; ++d)
        {
            if(created[i].first[d] > max[d])
            {
                max[d] = created[i].first[d];
                max_index[d] = i;
            }
            if(created[i].first[d] < min[d])
            {
                min[d] = created[i].first[d];
                min_index[d] = i;
            }
        }
    }
    //create initial polyhedron triangles from critical points MOOSE WARNING: verify that this is a valid convex hull--should always work for well-distributed pointsets, but need to prevent errors on weird input data
    Eigen::Vector3f &u = created[max_index[2]].first;
    Eigen::Vector3f &e = created[max_index[0]].first;
    Eigen::Vector3f &n = created[max_index[1]].first;
    Eigen::Vector3f &w = created[min_index[0]].first;
    Eigen::Vector3f &s = created[min_index[1]].first;
    Eigen::Vector3f &d = created[min_index[2]].first;
    Triangle *uen = new Triangle(u, e, n, 0);
    Triangle *unw = new Triangle(u, n, w, 0);
    Triangle *uws = new Triangle(u, w, s, 0);
    Triangle *use = new Triangle(u, s, e, 0);
    Triangle *des = new Triangle(d, e, s, 0);
    Triangle *dsw = new Triangle(d, s, w, 0);
    Triangle *dwn = new Triangle(d, w, n, 0);
    Triangle *dne = new Triangle(d, n, e, 0);
    uen->neighbors[0] = dne;    uen->neighbors[1] = unw;    uen->neighbors[2] = use;
    unw->neighbors[0] = dwn;    unw->neighbors[1] = uws;    unw->neighbors[2] = uen;
    uws->neighbors[0] = dsw;    uws->neighbors[1] = use;    uws->neighbors[2] = unw;
    use->neighbors[0] = des;    use->neighbors[1] = uen;    use->neighbors[2] = uws;
    des->neighbors[0] = use;    des->neighbors[1] = dsw;    des->neighbors[2] = dne;
    dsw->neighbors[0] = uws;    dsw->neighbors[1] = dwn;    dsw->neighbors[2] = des;
    dwn->neighbors[0] = unw;    dwn->neighbors[1] = dne;    dwn->neighbors[2] = dsw;
    dne->neighbors[0] = uen;    dne->neighbors[1] = des;    dne->neighbors[2] = dwn;
    //put our triangles into an array and connect them into a linked list
    Triangle *by_number[8] = { uen, unw, use, uws, dne, dwn, des, dsw };
    for(unsigned int i = 0; i < 7; ++i)
    {
        by_number[i]->next = by_number[i+1];
        by_number[i+1]->prev = by_number[i];
    }
    //sort critical point indices WARNING: max_index and min_index are invalidated by the sorting, but they aren't used again
    for(unsigned int n = 0; n < 5; ++n)
    {
        unsigned int least(n);
        for(unsigned int i = n + 1; i < 6; ++i)
            if(critical_index[i] < critical_index[least])
                least = i;
        if(least != n)
        {
            //swap critical_index[least] and critical_index[n]
            critical_index[least] ^= critical_index[n];
            critical_index[n] ^= critical_index[least];
            critical_index[least] ^= critical_index[n];
        }
    }
    //determine which triangle contains each point
    { // scope increase to prevent i from polluting enclosing scope
        unsigned int i = 0;
        for(unsigned int skip_index = 0; skip_index < 7; ++skip_index)
        {
            unsigned int last = (skip_index == 6 ? how_many : critical_index[skip_index]);
            for(; i < last; ++i)
            {
                //check the expected triangle
                if( by_number[created[i].second]->contains(created[i].first) )
                    ; //do nothing
                //check triangles sharing an edge with the expected triangle
                else if( by_number[(created[i].second ^ 1) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 1) & 7;
                else if( by_number[(created[i].second ^ 2) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 2) & 7;
                else if( by_number[(created[i].second ^ 4) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 4) & 7;
                //check triangles sharing a vertex with the expected triangle
                else if( by_number[(created[i].second ^ 6) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 6) & 7;
                else if( by_number[(created[i].second ^ 5) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 5) & 7;
                else if( by_number[(created[i].second ^ 3) & 7]->contains(created[i].first) )
                    created[i].second = (created[i].second ^ 3) & 7;
                //our point is contained in the triangle opposite to the expected triangle
                else
                    created[i].second = (~created[i].second) & 7;
                //increment the appropriate triangle's pointcount
                ++by_number[created[i].second]->pointcount;
            }
            ++i;
        }
    }
    //allocate space for points within triangles
    std::vector<std::pair<Eigen::Vector3f, unsigned int>*> dest_ptr; //dest_ptr[x] is the address where the next point inside destinations[x] should be copied
    for(unsigned int i = 0; i < 8; ++i)
    {
        dest_ptr.push_back(by_number[i]->points = (std::pair<Eigen::Vector3f, unsigned int>*) aligned_alloc( //align to cache so perform()'s multi-threaded RW operations can be optimized
            BowyerWatson::cache_line_size,
            ((by_number[i]->pointcount * BowyerWatson::point_size + BowyerWatson::cache_line_size - 1) / BowyerWatson::cache_line_size) * BowyerWatson::cache_line_size
        ));
    }
    //move points into triangles
    { // scope increase to prevent i from polluting enclosing scope
        unsigned int i = 0;
        for(unsigned int skip_index = 0; skip_index < 7; ++skip_index)
        {
            unsigned int last = (skip_index == 6 ? how_many : critical_index[skip_index]);
            for(; i < last; ++i)
                *(dest_ptr[created[i].second]++) = std::make_pair(created[i].first, (unsigned int)(-1));
            ++i;
        }
    }
    //return triangles
    return(std::make_pair(by_number[0], by_number[7]));
}

std::pair<BowyerWatson::Triangle*, BowyerWatson::Triangle*> BowyerWatson::create_cardinal_triangles(unsigned int how_many_minus_six) const
{
    std::vector<std::pair<Eigen::Vector3f, unsigned int> > created;
    generate_points(how_many_minus_six, created);
    return(create_cardinal_triangles(created));
}

std::pair<BowyerWatson::Triangle*, BowyerWatson::Triangle*> BowyerWatson::create_cardinal_triangles(std::vector<std::pair<Eigen::Vector3f, unsigned int> > &created) const
{
    unsigned int how_many = created.size();
    //create initial polyhedron
    Eigen::Vector3f u(0,0,1);
    Eigen::Vector3f e(1,0,0);
    Eigen::Vector3f n(0,1,0);
    Eigen::Vector3f w(-1,0,0);
    Eigen::Vector3f s(0,-1,0);
    Eigen::Vector3f d(0,0,-1);
    Triangle *uen = new Triangle(u, e, n, 0);
    Triangle *unw = new Triangle(u, n, w, 0);
    Triangle *uws = new Triangle(u, w, s, 0);
    Triangle *use = new Triangle(u, s, e, 0);
    Triangle *des = new Triangle(d, e, s, 0);
    Triangle *dsw = new Triangle(d, s, w, 0);
    Triangle *dwn = new Triangle(d, w, n, 0);
    Triangle *dne = new Triangle(d, n, e, 0);
    uen->neighbors[0] = dne;    uen->neighbors[1] = unw;    uen->neighbors[2] = use;
    unw->neighbors[0] = dwn;    unw->neighbors[1] = uws;    unw->neighbors[2] = uen;
    uws->neighbors[0] = dsw;    uws->neighbors[1] = use;    uws->neighbors[2] = unw;
    use->neighbors[0] = des;    use->neighbors[1] = uen;    use->neighbors[2] = uws;
    des->neighbors[0] = use;    des->neighbors[1] = dsw;    des->neighbors[2] = dne;
    dsw->neighbors[0] = uws;    dsw->neighbors[1] = dwn;    dsw->neighbors[2] = des;
    dwn->neighbors[0] = unw;    dwn->neighbors[1] = dne;    dwn->neighbors[2] = dsw;
    dne->neighbors[0] = uen;    dne->neighbors[1] = des;    dne->neighbors[2] = dwn;
    //create an array of our triangles and connect them into a linked list
    Triangle *by_number[8] = { uen, unw, use, uws, dne, dwn, des, dsw };
    for(unsigned int i = 0; i < 7; ++i)
    {
        by_number[i]->next = by_number[i+1];
        by_number[i+1]->prev = by_number[i];
    }
    //count the number of points per triangle
    for(unsigned int i = 0; i < how_many; ++i)
        ++by_number[created[i].second]->pointcount;
    //allocate space for points within triangles
    std::vector<std::pair<Eigen::Vector3f, unsigned int>*> dest_ptr; //dest_ptr[x] is the address where the next point inside destinations[x] should be copied
    for(unsigned int i = 0; i < 8; ++i)
    {
        dest_ptr.push_back(by_number[i]->points = (std::pair<Eigen::Vector3f, unsigned int>*) aligned_alloc(
            BowyerWatson::cache_line_size,
            ((by_number[i]->pointcount * BowyerWatson::point_size + BowyerWatson::cache_line_size - 1) / BowyerWatson::cache_line_size) * BowyerWatson::cache_line_size
        ));
    }
    //move points into triangles
    for(unsigned int i = 0; i < created.size(); ++i)
        *(dest_ptr[created[i].second]++) = std::make_pair(created[i].first, (unsigned int)(-1));
    //return triangles
    return(std::make_pair(by_number[0], by_number[7]));
}


BowyerWatson::Mesh BowyerWatson::get_mesh(Triangle *first, bool destructive, unsigned int total_points)
{
    //MOOSE WARNING: we don't use total_points, because this is an inefficient & crap algorithm that can't make use of it;
    //this function is to be rewritten later
    Mesh mesh;
    std::map<Eigen::Vector3f, std::vector<std::pair<unsigned int, unsigned int>>, PointComparator> points_to_triangles;
    unsigned int triangle_index = 0;
    for(Triangle *triangle = first; triangle != nullptr; triangle = triangle->next)
    {
        points_to_triangles[triangle->vertices[0]].push_back(std::make_pair(triangle_index, 0));
        points_to_triangles[triangle->vertices[1]].push_back(std::make_pair(triangle_index, 1));
        points_to_triangles[triangle->vertices[2]].push_back(std::make_pair(triangle_index, 2));
        ++triangle_index;
    }
    mesh.triangles.resize(triangle_index * 3);
    mesh.vertices.reserve(points_to_triangles.size() * 3);
    for(std::map<Eigen::Vector3f, std::vector<std::pair<unsigned int, unsigned int> >, PointComparator>::iterator iter = points_to_triangles.begin(); iter != points_to_triangles.end(); ++iter)
    {
        for(unsigned int i = 0; i < iter->second.size(); ++i)
            mesh.triangles[iter->second[i].first * 3 + iter->second[i].second] = mesh.vertices.size();
        mesh.vertices.push_back(iter->first[0]);
        mesh.vertices.push_back(iter->first[1]);
        mesh.vertices.push_back(iter->first[2]);
    }
    if(destructive)
    {
        Triangle *temp;
        while(first != 0)
        {
            temp = first->next;
            delete first;
            first = temp;
        }
    }
    return(mesh);
}




void BowyerWatson::relax_points(Triangle *start_triangle, std::vector<std::pair<Eigen::Vector3f, unsigned int> > &empty_vector, bool kill_triangles) const
{
    //calculate pre-centroids
    std::map<Eigen::Vector3f, std::pair<Eigen::Vector3f, float>, PointComparator> pre_centroids; //(point affected, (area-weighted sum of centroids, total area))
    for(Triangle *tri = start_triangle; tri != nullptr; tri = tri->next)
    {
        float area = (tri->vertices[1] - tri->vertices[0]).cross(tri->vertices[2] - tri->vertices[0]).norm() / 2.0f;
        Eigen::Vector3f centroid_times_area = area * (tri->vertices[0] + tri->vertices[1] + tri->vertices[2]) / 3.0f;
        for(unsigned int i = 0; i < 3; ++i)
        {
            std::pair<std::map<Eigen::Vector3f, std::pair<Eigen::Vector3f, float>, PointComparator>::iterator, bool> insert_data = pre_centroids.insert(
                std::make_pair(tri->vertices[i], std::make_pair(centroid_times_area, area))
            );
            if(!insert_data.second)
            {
                insert_data.first->second.first += centroid_times_area;
                insert_data.first->second.second += area;
            }
        }
    }
    //calculate centroids & place in empty_vector
    std::vector<std::pair<Eigen::Vector3f, unsigned int> > &created = empty_vector;
    created.clear();
    created.reserve(pre_centroids.size());
    for(std::map<Eigen::Vector3f, std::pair<Eigen::Vector3f, float>, PointComparator>::iterator cur_point = pre_centroids.begin(); cur_point != pre_centroids.end(); ++cur_point)
    {
        created.push_back(std::make_pair(cur_point->second.first / cur_point->second.second, 0));
        created.back().first /= created.back().first.norm();
        //guess which triangle contains the new point under the assumption that the critical points end up exactly at (1,0,0), (-1,0,0), (0,1,0), etc.
        //WARNING: the indices here refer to the by_number array of triangles, created later
        if(created.back().first[0] < 0) created.back().second |= 1;
        if(created.back().first[1] < 0) created.back().second |= 2;
        if(created.back().first[2] < 0) created.back().second |= 4;
    }
    //kill the triangles if requested
    if(kill_triangles)
        for(Triangle *tri = start_triangle; tri != nullptr; tri = tri->next)
            delete tri;
}





