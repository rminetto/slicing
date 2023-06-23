/*Authors: Rodrigo Minetto (UTFPR). */
/*         Jorge Stolfi (UNICAMP).  */
/*rodrigo.minetto@gmail.com*/

#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <set>
#include <unordered_map>
#include <limits>
#include <array>
#include <map>
#include <sys/stat.h>
#include <unistd.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>

#define GLM_FORCE_RADIANS

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/vector_angle.hpp"

#define GLM_FORCE_RADIANS
#define DEG_TO_RAD(x) (x*0.0174532925199f)
#define FILE_STL_BIN 0
#define FILE_STL_ASCII 1
#define FILE_AMF_ASCII 2

//#define DEBUG

using namespace std;

typedef CGAL::Quotient<CGAL::MP_Float> NT;

typedef CGAL::Cartesian<NT> Kernel;

typedef Kernel::Point_2 Point_2;

typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;

typedef Traits_2::Curve_2 Segment_2;

long intersections = 0;

double slicing_time = 0.0, loopclosure_time = 0.0;

/**************************************************************************
 *                        CLASS  point V3                                 * 
 **************************************************************************/
class v3 {

  public:
    
    v3 (float _x=0, float _y=0, float _z=0) : x(_x), y(_y), z(_z) {}

    float distTo (const v3 &pt) const { 
      return sqrt ( pow(fabs(x-pt.x), 2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z), 2.0) ); 
    }

    array<float,3> getCoords() { 
      array<float,3> c = {{x, y, z}}; 
      return c; 
    }

    float dotproduct (const v3 &v) const { 
      return (x*v.x + y*v.y + z*v.z); 
    }

    void transform (const glm::mat4 &mat) { 
      glm::vec4 v = glm::vec4(x, y, z, 1.0);
      glm::vec4 vt = mat*v; 
      x = (vt.x); y = (vt.y); z = (vt.z);
    }
    
    v3& operator-=(const v3 &pt) { 
      x = (x-pt.x); 
      y = (y-pt.y); 
      z = (z-pt.z); 
      return *this;    
    }

    v3 operator-(const v3 &pt) { 
      return v3 ((x-pt.x), (y-pt.y), (z-pt.z)); 
    }
    
    v3 operator+(const v3 &pt) { 
      return v3 ((x+pt.x), (y+pt.y), (z+pt.z)); 
    }
    
    v3 operator/(float a) { 
      return v3 ((x/a), (y/a), (z/a)); 
    }
    
    v3 operator*(float a) { 
      return v3 ((x*a), (y*a), (z*a)); 
    }
    
    bool operator<(const v3 &pt) const { 
      return z < pt.z; 
    }
    
    bool operator>(const v3 &pt) const { 
      return z > pt.z; 
    }
    
    bool operator==(const v3 &pt) const {
      return distTo(pt) < 0.005; 
    }
    
    bool operator!=(const v3 &pt) const { 
      return distTo(pt) > 0.005; 
    }
    
    float normalize() const { 
      return sqrt(x*x+y*y+z*z); 
    }

    string getLabel() const {
      stringstream ss;
      ss << x << "|" << y << "|" << z;
      return ss.str();
    }
    
    friend ostream& operator<<(ostream& os, const v3& v) {
      os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
      return os;
    }

  public:

    float x;
    float y;
    float z;
};

v3 operator-(const v3 &a, const v3 &b) {return v3((a.x-b.x), (a.y-b.y), (a.z-b.z)); }

v3 operator+(const v3 &a, const v3 &b) {return v3((a.x+b.x), (a.y+b.y), (a.z+b.z)); }

/**************************************************************************
 *                        CLASS  LineSegment                              * 
 **************************************************************************/
class LineSegment {

  public: 

    LineSegment (v3 p0=v3(), v3 p1=v3(), int i=0) { 
      v[0] = p0; 
      v[1] = p1;
      index = i;
      vertical = false;  
      if ((v[1].x - v[0].x) != 0) {
        a = (v[1].y - v[0].y)/(v[1].x - v[0].x);
        b = (v[0].y - (a * v[0].x));
      }
      else {
        vertical = true;
      }
    }

    bool operator==(const LineSegment &ls) const { 
      return ((v[0] == ls.v[0]) && (v[1] == ls.v[1])); 
    }

    friend ostream& operator<<(ostream& os, const LineSegment& ls) {
      os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")";
      return os;
    }

  public:

    v3 v[2];
    double a;
    double b;
    bool vertical;
    int index;
}; 


/**************************************************************************
 *                                  HASH                                  * 
 **************************************************************************/
template<typename T> inline void hash_combine(size_t &seed, const T &v) {
  hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct HashV3 {
  size_t operator() (const v3 &v) const {
    size_t h = hash<float>()(v.x);
    hash_combine (h, v.y);
    hash_combine (h, v.z);
    return h;
  }
};

typedef unordered_map<v3, vector<v3>, HashV3> PointMesh;


/**************************************************************************
 *                        CLASS  Triangle                                 * 
 **************************************************************************/
class Triangle    {

  public:
  
    Triangle(v3 n, v3 v0, v3 v1, v3 v2) : normal(n) {
      v[0] = v0; 
      v[1] = v1; 
      v[2] = v2;
      zMin = +99999999.9; 
      zMax = -99999999.9;
      setZMin(v0.z); setZMin(v1.z); setZMin(v2.z);
      setZMax(v0.z); setZMax(v1.z); setZMax(v2.z);
    }

    void setZMin (float z) { 
      if (z < zMin) {
        zMin = z; 
      }
    }

    void setZMax (float z) { 
      if (z > zMax) {
        zMax = z; 
      }
    }

    Triangle& operator-=(const v3 &pt) { 
      v[0] -= pt; 
      v[1] -= pt;  
      v[2] -= pt; 
      return *this;
    }

    bool operator<(const Triangle &t) { 
       return zMin < t.zMin; 
    }

    friend ostream& operator<<(ostream& os, const Triangle& t) {
      os << "V0: (" << t.v[0] << "); V1: (" << t.v[1] << "); V2: (" << t.v[2] << ")";
      return os;
    }
  
  public:

    v3 v[3];
    v3 normal;
    float zMin;
    float zMax;
};  

/**************************************************************************
 *                      CLASS  TriangleMesh                               * 
 **************************************************************************/
class TriangleMesh {

  public:
  
    TriangleMesh() : bottomLeftVertex(999999,999999,999999), upperRightVertex(-999999,-999999,-999999) { meshSize = 0;}

    size_t size() const {
      return meshSize;
    }

    void push_back(Triangle &t) {
      meshSize++;
      vTriangle.push_back(t);
      for (size_t i = 0; i < 3; ++i) {
        if (t.v[i].x < bottomLeftVertex.x) { bottomLeftVertex.x = t.v[i].x; }
        if (t.v[i].y < bottomLeftVertex.y) { bottomLeftVertex.y = t.v[i].y; }
        if (t.v[i].z < bottomLeftVertex.z) { bottomLeftVertex.z = t.v[i].z; }
        if (t.v[i].x > upperRightVertex.x) { upperRightVertex.x = t.v[i].x; }
        if (t.v[i].y > upperRightVertex.y) { upperRightVertex.y = t.v[i].y; }
        if (t.v[i].z > upperRightVertex.z) { upperRightVertex.z = t.v[i].z; }
      }
    }

    v3 meshAABBSize() const {
        return v3 ( upperRightVertex.x - bottomLeftVertex.x, 
                    upperRightVertex.y - bottomLeftVertex.y, 
                    upperRightVertex.z - bottomLeftVertex.z );
    }

    const vector<Triangle>& getvTriangle() const { return vTriangle; }

    v3 getBottomLeftVertex() const { return bottomLeftVertex; }

    v3 getUpperRightVertex() const { return upperRightVertex; }

  public:

    int meshSize;
    v3 bottomLeftVertex;
    v3 upperRightVertex;
    vector<Triangle> vTriangle;
};

/*Functions prototypes: */
void export_svg_no_chaining (string fileName, vector<LineSegment> Segments[], int k, const v3 &aabbSize);

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*******************************************Contour orientation***************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

typedef struct _rgb {
    int r, g, b;
} rgb;

typedef struct _bounding_box {
  double xMin;
  double xMax;
  double yMin;
  double yMax;
  bool firstPoint;
} bounding_box;

typedef struct _contour {
   bool external;
   bool clockwise;
   vector<v3> P;
} contour;

/*-----------------------------------------------------------------------*/
void record_polygons (vector<contour> polygons, FILE *file) {
  for (int i = 0; i < polygons.size(); i++) {
    vector<v3> Pi = polygons.at(i).P;
    int size = (int)(Pi.size()) - 1;
    fprintf(file, "%d\n", size);
    for (int p = 0; p < size; p++) {
       v3 v = Pi.at(p);
       fprintf(file, "%f %f\n", v.x, v.y);
    }
  }
}

/*-----------------------------------------------------------------------*/
bool is_inside (LineSegment line, v3 point) {
  double maxX = (line.v[0].x > line.v[1].x) ? line.v[0].x : line.v[1].x;
  double minX = (line.v[0].x < line.v[1].x) ? line.v[0].x : line.v[1].x;
  double maxY = (line.v[0].y > line.v[1].y) ? line.v[0].y : line.v[1].y;
  double minY = (line.v[0].y < line.v[1].y) ? line.v[0].y : line.v[1].y;
  if ((point.x >= minX && point.x <= maxX) && (point.y >= minY && point.y <= maxY)) {
    return true;
  }
  return false;
}

/*-----------------------------------------------------------------------*/
bool ray_intersect (LineSegment ray, LineSegment side) {
  v3 intersectPoint;
  /* If both vectors aren't from the kind of x=1 lines then go into: */
  if (!ray.vertical && !side.vertical) {
    /* Check if both vectors are parallel. If they are parallel then no intersection point will exist: */
    if (ray.a - side.a == 0) {
      return false;
    }
    intersectPoint.x = ((side.b - ray.b) / (ray.a - side.a)); 
    intersectPoint.y = side.a * intersectPoint.x + side.b; 
  }
  else if (ray.vertical && !side.vertical) {
    intersectPoint.x = ray.v[0].x;
    intersectPoint.y = side.a * intersectPoint.x + side.b;
  }
  else if (!ray.vertical && side.vertical) {
    intersectPoint.x = side.v[0].x;
    intersectPoint.y = ray.a * intersectPoint.x + ray.b;
  }
  else {
    return false;
  }
  if (is_inside(side, intersectPoint) && is_inside(ray, intersectPoint)) {
    return true;
  }
  return false;
}

/*-----------------------------------------------------------------------*/
LineSegment create_ray (v3 point, bounding_box bb, int index) {
  /* Create outside point: */
  double epsilon = (bb.xMax - bb.xMin) / 100.0;
  v3 outside (bb.xMin - epsilon, bb.yMin);
  LineSegment v (outside, point, index);
  return v;
}

/*-----------------------------------------------------------------------*/
bounding_box create_bounding_box () {
  bounding_box bb;
  bb.xMax = std::numeric_limits<double>::min();
  bb.xMin = std::numeric_limits<double>::max();
  bb.yMax = std::numeric_limits<double>::min();
  bb.yMin = std::numeric_limits<double>::max();
  return bb;
}

/*-----------------------------------------------------------------------*/
void update_bounding_box (v3 point, bounding_box *bb) {
  /* Setting the bounding box: */
  if (point.x > bb->xMax) {
    bb->xMax = point.x;
  }
  else if (point.x < bb->xMin) {
    bb->xMin = point.x;
  }
  if (point.y > bb->yMax) {
    bb->yMax = point.y;
  }
  else if (point.y < bb->yMin) {
    bb->yMin = point.y;
  }
}

/*-----------------------------------------------------------------------*/
bool insided_bounding_box (v3 point, bounding_box bb) {
  if ( (point.x < bb.xMin) || (point.x > bb.xMax) || (point.y < bb.yMin) || (point.y > bb.yMax) ) {
     return false;
  }
  return true;
}

/*-----------------------------------------------------------------------*/
bool contains (v3 point, bounding_box bb, vector<LineSegment> sides, int index) {
  if (insided_bounding_box(point, bb)) {
    LineSegment ray = create_ray (point, bb, index);
    int intersection = 0;
    for (int i = 0; i < sides.size(); i++) {
      if ((sides.at(i).index != index) && ray_intersect(ray, sides.at(i))) {
         intersection++;
      }
    }
    /* If the number of intersections is odd, then the point is inside the polygon: */
    if ((intersection % 2) == 1) {
      return true;
    }
  }
  return false;
}

/*-----------------------------------------------------------------------*/
void add_point (v3 p1, v3 p2, vector<LineSegment> &t, bounding_box *bb, bool first, int index) {
   if (first) {
     update_bounding_box (p1, bb);
   }
   update_bounding_box (p2, bb);
   LineSegment line (p1, p2, index);
   t.push_back(line);
}

/*-----------------------------------------------------------------------*/
/*Function to (re)orient the contour clockwise and counter-clockwise.*/
void ray_casting (vector<contour> &polygons) {

  vector<LineSegment> segments;
  
  bounding_box bb = create_bounding_box ();

  /*Creating the line segments of each contour: */
  for (int i = 0; i < polygons.size(); i++) {
    double area = 0.0;
    vector<v3> Pi = polygons.at(i).P;
    for (int j = 1; j < Pi.size(); j++) {
      v3 p0 = Pi.at(j-1);
      v3 p1 = Pi.at(j+0);
      area += (p0.x * p1.y - p0.y * p1.x);
      add_point (p0, p1, segments, &bb, (j == 1 ? true : false), i);
      if (j == Pi.size()-1) {
         add_point (p1, Pi.at(0), segments, &bb, (j == 1 ? true : false), i);
         area += (p1.x * Pi.at(0).y - p1.y * Pi.at(0).x);
      }
    }
    area /= 2.0;
    if (area < 0.0) {
      polygons.at(i).clockwise = true;
    }
    else {
      polygons.at(i).clockwise = false;
    }
  }

  /*Using the point in polygon algorithm to test the first segment of each contour: */
  for (int i = 0; i < polygons.size(); i++) {
    vector<v3> Pi = polygons.at(i).P;
    if (contains (Pi.at(0), bb, segments, i)) {
      /*Internal contour: */
      polygons.at(i).external = false;
    }
    else {
      /*External contour: */
      polygons.at(i).external = true;
    }

    /*Reversing contours: */
    if (polygons.at(i).external && polygons.at(i).clockwise) {
       std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
       polygons.at(i).clockwise = false;
    }
    else if (!polygons.at(i).external && !polygons.at(i).clockwise) {
       std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
       polygons.at(i).clockwise = true;
    }
  }
  segments.clear();
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/**********************************************Incremental********************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

/* Assumes that {P} is a list of {k} strictly increasing {Z} coordinates. 
   Returns an integer {p} such that {P[p-1] < zMin < P[p]}. As special cases, 
   if {zMin < P[0]} returns 0; if {zMin > P[k-1]} returns {k}. */
int IncrementalSlicing_binary_search (float zMin, vector<float> P) {
  int k = P.size();
  assert(k >= 1);
  if (zMin >= P[k-1]) { return k; }
  /* Binary search: */
  int l = -1; /* Inferior Z index. */
  int r = k;  /* Superior Z index. */
  while (r - l > 1) {
    /* At this point, {zMin} is between {P[l]} and {P[r]}. */
    int m = (l + r)/2;
    assert((0 <= m) && (m < k));
    if (zMin >= P[m]) {
      l = m;
    } else {
      r = m;
    }
  }
  return r;
}

/*-----------------------------------------------------------------------*/
/*Structures*/
typedef struct node {
  Triangle t;
  struct node *next;
  struct node *prev;
} Mesh_Triangle_Node_t;

typedef struct _list {
   Mesh_Triangle_Node_t *head;
   Mesh_Triangle_Node_t *tail;
} Mesh_Triangle_List_t;

/*-----------------------------------------------------------------------*/
Mesh_Triangle_List_t* Mesh_Triangle_List_create (void) {
  Mesh_Triangle_List_t *L = (Mesh_Triangle_List_t *)malloc(sizeof(Mesh_Triangle_List_t));
  L->head = NULL;
  L->tail = NULL;
  return L;
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_insert (Triangle t, Mesh_Triangle_List_t *L) {
  Mesh_Triangle_Node_t *node = (Mesh_Triangle_Node_t *)malloc(sizeof(Mesh_Triangle_Node_t));
  node->t = t;
  node->next = L->head;
  node->prev = NULL;
  if (L->head == NULL) {
    /*New head*/
    L->head = L->tail = node;
  }
  else {
    L->head->prev = node;
    L->head = node;
  }
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2) {
  if ( (L1->head != NULL) && (L2->head != NULL) ) {
    L1->tail->next = L2->head;
    L2->head->prev = L1->tail;
    L1->tail = L2->tail;;
  }
  else if (L2->head != NULL) {
    L1->head = L2->head;
    L1->tail = L2->tail;
  }
}

/*-----------------------------------------------------------------------*/
Mesh_Triangle_Node_t* Mesh_Triangle_List_remove (Mesh_Triangle_List_t *L, Mesh_Triangle_Node_t *node) {
  if ((node->prev == NULL) && (node->next == NULL)) {
    free (node);
    L->head = NULL;
    L->tail = NULL;
    return NULL;
  }
  else if (node->prev == NULL) {
    node->next->prev = NULL;
    L->head = node->next;
    free (node);
    return L->head;
  }
  else if (node->next == NULL) {
    node->prev->next = NULL;
    L->tail = node->prev;
    free (node);
    return NULL;
  }
  else {
    Mesh_Triangle_Node_t *next = node->next;
    node->next->prev = node->prev;
    node->prev->next = next;
    free (node);
    return next;
  }
}

/*-----------------------------------------------------------------------*/
v3 R3_Mesh_Side_slice (v3 vi, v3 vj, float Z) {
   double dx = vj.x - vi.x;
   double dy = vj.y - vi.y;
   double dz = vj.z - vi.z;
   assert(dz != 0);
   double frac = (Z - vi.z)/dz;
   float xint = (float)(frac*dx + (double)vi.x);
   float yint = (float)(frac*dy + (double)vi.y);
   //return (v3){ .x = xint, .y = yint, .z = Z };
   return (v3){xint, yint, Z};
}

/*-----------------------------------------------------------------------*/
LineSegment R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z) {
   assert((t->t.zMin < Z) && (t->t.zMax > Z));
   int np = 0; /* Number of segment endpoints found */
   LineSegment seg;
   for (int i = 0; i < 3; i++) {
      /* Get side {i} of triangle: */
      int j = (i == 2 ? 0 : i+1);
      v3 vi = (t->t.v[i]);
      v3 vj = (t->t.v[j]);
      /* Check for intersection of plane with {vi--vj}. */
      /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
      float vzMin = (vi.z < vj.z ? vi.z : vj.z);
      float vzMax = (vi.z > vj.z ? vi.z : vj.z);
      if ((vzMin <= Z) && (vzMax > Z)) {
         v3 p = R3_Mesh_Side_slice (vi, vj, Z);
         assert(np < 2);
         seg.v[np] = p;
         np++;
      }
   }
   assert(np == 2);
   return seg; 
}

/*-----------------------------------------------------------------------*/
LineSegment R3_Mesh_Triangle_slice (Triangle t, float Z) {
   assert((t.zMin < Z) && (t.zMax > Z));
   int np = 0; /* Number of segment endpoints found */
   LineSegment seg;
   for (int i = 0; i < 3; i++) {
      /* Get side {i} of triangle: */
      int j = (i == 2 ? 0 : i+1);
      v3 vi = (t.v[i]);
      v3 vj = (t.v[j]);
      /* Check for intersection of plane with {vi--vj}. */
      /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
      float vzMin = (vi.z < vj.z ? vi.z : vj.z);
      float vzMax = (vi.z > vj.z ? vi.z : vj.z);
      if ((vzMin <= Z) && (vzMax > Z)) {
         v3 p = R3_Mesh_Side_slice (vi, vj, Z);
         assert(np < 2);
         seg.v[np] = p;
         np++;
      }
   }
   assert(np == 2);
   return seg; 
}

/*----------------------------------------------------------------------*/
/* Assumes that {P[0..k-1]} is a list of {k} strictly increasing {Z} 
  coordinates. Returns a vector of {k+1} lists {L[0..k]} such that {L[p]} 
  contains all triangles of the {mesh} that have {zMin} between {P[p-1]} 
  and {P[p]}, assuming that {P[-1] = -oo} and {P[k] = +oo}. If {delta > 0}, 
  assumes that {P[p]-P[p-1] = delta} for {p} in {1..k-1}. If {srt} is true, 
  assumes that the triangles are already sorted by increasing {zMin}. */
Mesh_Triangle_List_t** IncrementalSlicing_buildLists (bool srt, double delta, const TriangleMesh *mesh, vector<float> P) {

  int k = P.size(); /* Number of planes. */

  Mesh_Triangle_List_t **L = (Mesh_Triangle_List_t **)malloc((k+1) * sizeof(Mesh_Triangle_List_t *));

  for (size_t p = 0; p <= k; p++) { L[p] = Mesh_Triangle_List_create(); }

  const vector<Triangle> &T = mesh->getvTriangle();

  int n = T.size(); /* Number of triangles. */

  if (delta > 0.0) {
    /* Uniform slicing - compute list index: */
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      Triangle t = *it;
      int p;
      if (t.zMin < P[0]) {
        p = 0;
      }
      else if (t.zMin > P[k-1]) {
        p = k;
      }
      else {
        p = floor((t.zMin - P[0])/delta) + 1;
      }
      Mesh_Triangle_List_insert (t, L[p]);
    }
  } else if (srt) {
    /* Slicing of a pre-sorted mesh - merge {zMin}s and {P}: */
    auto it = T.begin();
    auto itEnd = T.end();
    double zprev = -INFINITY;
    for (int p = 0; p <= k; p++) {
      float Zp = (p < k ? P[k] : +INFINITY);
      const Triangle t = *it;
      assert(t.zMin >= zprev);
      while ((it != itEnd) && (t.zMin < Zp)) {
        Mesh_Triangle_List_insert (t, L[p]);
        zprev = t.zMin;
        it++;
      }
    }
  } else {
    /* General case: */
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      const Triangle t = *it;
      int p = IncrementalSlicing_binary_search(t.zMin, P);
      assert((p >= 0) && (p <= k));
      Mesh_Triangle_List_insert (t, L[p]);
    }
  }
  return L;
}

/*-----------------------------------------------------------------------*/
/*Gets an arbitrary segment from {H}, removes it from {H} and returns it as a trivial chain. */
vector<v3> IncrementalStartLoop(vector<PointMesh> &H) {
   vector<v3> P;
   auto it = (H[0]).begin();
   v3 u = (*it).first;
   vector<v3> vw = (*it).second;
   v3 v = vw.at(0);
   P.push_back(u);
   P.push_back(v);
   (H[0][u]).erase(std::remove((H[0][u]).begin(), (H[0][u]).end(), v), (H[0][u]).end());
   if (H[0][u].size() == 0) { (H[0]).erase(u); }
   (H[0][v]).erase(std::remove((H[0][v]).begin(), (H[0][v]).end(), u), (H[0][v]).end());
   if (H[0][v].size() == 0) { (H[0]).erase(v); }
   return P;
}

/*-----------------------------------------------------------------------*/
/*Extends the chain {P} wih segments from {H}, removing them, while possible. */
void IncrementalExtendLoop(vector<v3> &P, vector<PointMesh> &H) { 
  int index = 0;
  int n = P.size();
  v3 first = P.front();  
  v3 current = P.back(); 
  v3 last;
         
  /* Collect other vertices: */
  while (true) {
    auto it = (H[0]).find(current);
    if (it == (H[0]).end()) { /*Vertex {current} is a dead end:*/ break; }
    v3 key1 = (*it).first; assert(key1 == current);  /*Paranoia check.*/
            
    /*Get {next}, the first unused neighbor of {current}:*/
    vector<v3> vw = (*it).second; /*Unused neighbors of {current}.*/
    assert (vw.size() != 0); 
    v3 next = vw.at(0); /*First unused neighbor of {current}.*/

    /*Append the segment {(current,next)} to {P} and delete from {H}:*/
    P.push_back(next);

    /*Remove the segment {(current,next)} from {H}:*/
    (H[0][current]).erase(std::remove((H[0][current]).begin(), (H[0][current]).end(), next), (H[0][current]).end());
    if (H[0][current].size() == 0) { (H[0]).erase(current); } 
    (H[0][next]).erase(std::remove((H[0][next]).begin(), (H[0][next]).end(), current), (H[0][next]).end());
    if (H[0][next].size() == 0) { (H[0]).erase(next); } 

    if (next == first) { /*We closed a loop:*/ break; }

    /*Move on:*/
    current = next;
  }
}

/*Reverses the chain {P}.*/
void IncrementalReverseLoop(vector<v3> &P) { 
  std::reverse(P.begin(),P.end());
}

/*-----------------------------------------------------------------------
  Groups the segments {segs} into zero or more open or closed polygonal chains.
A chain with {n} segments has {n+1} vertices (of type {v3}). There may be
repeated vertices.  The chain is closed iff the last vertex is equal to the
first. Vertices have their {x} and {y} coordinates rounded to a multiple of a
certain {eps}. The output will have one open chain for every pair of vertices
of odd degree. The rest of the segments will be grouped into closed chains.
All chains are returned in the {Polygons} list. */
void ContourConstruction (vector<LineSegment> segs, vector<contour> polygons[], int plane) {

   bool verbose = false;

   clock_t contour_begin = clock();

   /*Creating the hash table.*/
   vector<PointMesh> H(1);

   /*Rounding vertices and filling the hash table.*/
   double eps = 1/128.0;
   for (std::vector<LineSegment>::iterator i = segs.begin(); i != segs.end(); i++) {
      LineSegment q = *i;
      q.v[0].x = round(q.v[0].x / eps) * eps;
      q.v[0].y = round(q.v[0].y / eps) * eps;
      q.v[0].z = plane;
      q.v[1].x = round(q.v[1].x / eps) * eps;
      q.v[1].y = round(q.v[1].y / eps) * eps;
      q.v[1].z = plane;
      if (q.v[0].distTo(q.v[1]) > 0.0001) {
         (H[0][q.v[0]]).push_back(q.v[1]);
         (H[0][q.v[1]]).push_back(q.v[0]);
      }
   }

   /* Count vertices by degree: */
   if (verbose) {
     int degmax = 10;
     int ctdeg[degmax+1];
     for (int deg = 0; deg <= degmax; deg++) { ctdeg[deg] = 0; }
     for (auto i = (H[0]).begin(); i != (H[0]).end(); i++) {
        vector<v3> L = (*i).second;
        int deg = L.size();
        if (deg > degmax) { deg = degmax; }
        ctdeg[deg]++;
     }
     assert(ctdeg[0] == 0);
     bool closedSlice = true;
     for (int deg = 1; deg <= degmax; deg++) { 
       if (((deg % 2) != 0) && (ctdeg[deg] > 0)) { closedSlice = false; }
       if ((verbose || (deg != 2)) && (ctdeg[deg] != 0))
         { cout << "there are " << ctdeg[deg] << " vertices of degree " << deg << " on plane " << plane << endl; }
     }
     if (!closedSlice) { cout << "** contours of plane " << plane << " are not closed" << endl; }
   }

   /*Contour construction.*/
   bool maximal = true;
   while (!(H[0]).empty()) {
     if (maximal) {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       if (P.front() != P.back()) { //Chain {P} is open
         IncrementalReverseLoop(P);
         IncrementalExtendLoop(P,H);
       }
       polygons[plane].push_back({false, false, P});
     }
     else {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       polygons[plane].push_back({false, false, P});
     }
   }
   clock_t contour_end = clock();
   loopclosure_time += double(contour_end - contour_begin)/CLOCKS_PER_SEC;
}

/*-----------------------------------------------------------------------*/
void IncrementalSlicing (const TriangleMesh *mesh, vector<float> P, float delta, bool srt, vector<contour> polygons[], bool chaining, bool orienting) {

  /*Slicing*/
  clock_t slice_begin = clock();

  int k = P.size();

  vector<LineSegment> segs[k];

  /* Classify triangles by the plane gaps that contain their {zMin}: */
  Mesh_Triangle_List_t **L = IncrementalSlicing_buildLists (srt, delta, mesh, P);
  /* Now perform a plane sweep from bottom to top: */

  Mesh_Triangle_List_t *A = Mesh_Triangle_List_create(); /* Active triangle list. */
  for (int p = 0; p < k; p++) {
    /* Add triangles that start between {P[p-1]} and {P[p]}: */
    Mesh_Triangle_List_union (A, L[p]);
    /* Scan the active triangles: */
    Mesh_Triangle_Node_t *aux = A->head;
    while (aux != NULL) {
      Mesh_Triangle_Node_t *next = aux->next;
      if (aux->t.zMax < P[p]) {
        /* Triangle is exhausted: */
        Mesh_Triangle_List_remove (A, aux);
      } else {
        /* Compute intersection: */
        if ((aux->t.zMin < P[p]) && (aux->t.zMax > P[p])) {
          LineSegment seg = R3_Mesh_Triangle_slice (aux, P[p]);
          segs[p].push_back(seg);
          intersections++;
        }
      }
      aux = next;
    }
  }
  free(L);
  clock_t slice_end = clock();
  slicing_time = double(slice_end - slice_begin)/CLOCKS_PER_SEC;
  /*End-Slicing*/

  if (chaining) {
    /*Contour construction:*/
    for (size_t p = 0; p < k; p++) {
      if (!segs[p].empty()) {
         ContourConstruction (segs[p], polygons, p);
         #ifdef DEBUG
           char fname[256];
           sprintf(fname, "slice_%03d.txt", (int)p);
           FILE *fslice = fopen(fname, "w");
           fprintf(fslice, "----------- Segmentos ---------------\n");
           for (int ii = 0; ii < segs[p].size(); ii++) {
              LineSegment ss = segs[p].at(ii); 
              fprintf(fslice, "%f %f   %f %f\n", ss.v[0].x, ss.v[0].y, ss.v[1].x, ss.v[1].y);           
           }
           fprintf(fslice, "---------- End Segmentos --------------\n");
           record_polygons (polygons[p], fslice);
           fclose(fslice);
         #endif
         if (orienting) {   
           ray_casting (polygons[p]);
         }
         segs[p].clear(); 
      }
    }
    /*End construction.*/
  }
  else {
    export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
  }
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/**************************************************Trivial********************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
v3 TrivialClosure_find_B (vector<LineSegment>& segs, v3 u, int plane) {
  double dist = +INFINITY;
  v3 sp;
  int piter = 1, pitem = -1;
  for (vector<LineSegment>::iterator j = segs.begin(); j != segs.end(); j++) {
    v3 p0 = (*j).v[0];
    v3 p1 = (*j).v[1];
    if ((p0.x == u.x) && (p0.y == u.y)) {
      segs.erase(j);
      return p1;
    } else if ((p1.x == u.x) && (p1.y == u.y)) {
      segs.erase(j);
      return p0;
    } else {
      double d0 = sqrt ((p0.x - u.x) * (p0.x - u.x) + (p0.y - u.y) * (p0.y - u.y));
      double d1 = sqrt ((p1.x - u.x) * (p1.x - u.x) + (p1.y - u.y) * (p1.y - u.y));
      if (d0 < dist) { dist = d0; sp = p0; pitem = piter; }
      if (d1 < dist) { dist = d1; sp = p1; pitem = piter; }
    }
    piter++;
  }
  if ((pitem != -1) && (dist <= 1)) {
    segs.erase(segs.begin() + pitem);
    return sp;
  }
  else {
    cout << "Distance: " << dist << ", plane : "<< plane << endl; 
    return {+INFINITY,+INFINITY,+INFINITY};
  }
}

/*-----------------------------------------------------------------------*/
v3 TrivialClosure_find_A (vector<LineSegment>& segs, v3 u, int plane) {
  for (vector<LineSegment>::iterator j = segs.begin(); j != segs.end(); ++j) {
    v3 p0 = (*j).v[0];
    v3 p1 = (*j).v[1];
    if ((p0.x == u.x) && (p0.y == u.y)) {
      segs.erase(j);
      return p1;
    } else if ((p1.x == u.x) && (p1.y == u.y)) {
      segs.erase(j);
      return p0;
    }
  }
  return {+INFINITY,+INFINITY,+INFINITY};
}

/*-----------------------------------------------------------------------*/
void TrivialLoopClosure (vector<LineSegment> segs, vector<contour> polygons[], int plane) {

  while (segs.size() > 0) {

    /* Get another contour: */
    vector<v3> P;

    /* Get the first segment: */
    v3 last;
    v3 current;
    v3 prev;
    {
      std::vector<LineSegment>::iterator i = segs.begin();
      last = (*i).v[0];
      current = (*i).v[1];
      {
         P.push_back (last);
      }
      prev = last;
      segs.erase(i);
    }

    /* Get additional segments until loop is closed: */
    bool open = false;
    do {
      /* Find and delete another segment with endpoint {current}, advance {prev,current} */
      v3 next = TrivialClosure_find_A (segs, current, plane);
      //v3 next = TrivialClosure_find_B (segs, current, plane);
      if ((next.x == +INFINITY) || (next.y == +INFINITY)) {
        /* Open contour?! */
        open = true;
        break;
      }
      {
         P.push_back (current);
      }
      prev = current;
      current = next;
    } while ((current.x != last.x) || (current.y != last.y));
    if (!open) {
      P.push_back(last);
    } 
    polygons[plane].push_back({false, false, P});
  }
}

/*-----------------------------------------------------------------------*/
void TrivialSlicing (const TriangleMesh *mesh, vector<float> P, vector<contour> polygons[], bool chaining, bool orienting) {

  /*Slicing*/
  clock_t slice_begin = clock();

  int k = P.size(); /* Number of planes. */

  const vector<Triangle> &T = mesh->getvTriangle();

  vector<LineSegment> segs[k];

  /* Enumerate the slicing planes: */
  for (int p = 0; p < k; p++) {
    /* Enumerate all triangles of the mesh:*/
    int intersections_per_plane = 0;
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      Triangle t = *it; /*Current triangle.*/
      /*Does the plane intersect the triangle?:*/
      if ((t.zMin < P[p]) && (t.zMax > P[p])) {
        /* Compute and save the intersection: */
        LineSegment seg = R3_Mesh_Triangle_slice (t, P[p]);
        seg.v[0].x = round(seg.v[0].x * 100.0) / 100.0;
        seg.v[0].y = round(seg.v[0].y * 100.0) / 100.0;
        seg.v[1].x = round(seg.v[1].x * 100.0) / 100.0;
        seg.v[1].y = round(seg.v[1].y * 100.0) / 100.0;
        if (seg.v[0].distTo(seg.v[1]) > 0.0001) {
           segs[p].push_back(seg);
        }
        intersections++;
        intersections_per_plane++;
      }
    }
  }
  clock_t slice_end = clock();
  slicing_time = double(slice_end - slice_begin)/CLOCKS_PER_SEC;
  /*End-Slicing*/
       
  if (chaining) {
    /*Loop-Closure:*/
    clock_t contour_begin = clock();
    for (size_t p = 0; p < k; p++) {
      if (!segs[p].empty()) {
        TrivialLoopClosure (segs[p], polygons, p);
        if (orienting) {
           ray_casting (polygons[p]);
        }
        #ifdef DEBUG
           char fname[256];
           sprintf(fname, "slice_%03d.txt", (int)p);
           FILE *fslice = fopen(fname, "w");
           record_polygons (polygons[p], fslice);
           fclose(fslice);
         #endif
        segs[p].clear(); 
      }
    }
    clock_t contour_end = clock();
    loopclosure_time = double(contour_end - contour_begin)/CLOCKS_PER_SEC;
    /*End-Loop-Closure*/
  }
  else {
    export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
  }
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/***************************************************Park**********************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

/*-----------------------------------------------------------------------*/
void ParkLoopClosure (vector<Segment_2> in, vector<contour> polygons[]) {
   std::list<Point_2> pts;
   CGAL::compute_intersection_points (in.begin(), in.end(), std::back_inserter(pts), true);
   in.clear();
   pts.clear();
}

/*-----------------------------------------------------------------------*/
void ParkSlicing (const TriangleMesh *mesh, vector<float> P, vector<contour> polygons[], bool chaining) {

  /*Slicing*/
  clock_t slice_begin = clock();

  int k = P.size(); /* Number of planes. */

  const vector<Triangle> &T = mesh->getvTriangle();

  vector<LineSegment> segs[k];

  /* Enumerate the slicing planes: */
  for (int p = 0; p < k; p++) {
    /* Enumerate all triangles of the mesh:*/
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      Triangle t = *it; /*Current triangle.*/
      /*Does the plane intersect the triangle?:*/
      if ((t.zMin < P[p]) && (t.zMax > P[p])) {
        /* Compute and save the intersection: */
        LineSegment seg = R3_Mesh_Triangle_slice (t, P[p]);
        seg.v[0].x = round(seg.v[0].x * 100.0) / 100.0;
        seg.v[0].y = round(seg.v[0].y * 100.0) / 100.0;
        seg.v[1].x = round(seg.v[1].x * 100.0) / 100.0;
        seg.v[1].y = round(seg.v[1].y * 100.0) / 100.0;
        if (seg.v[0].distTo(seg.v[1]) > 0.0001) {
           segs[p].push_back(seg);
        }
        intersections++;
      }
    }
  }
  clock_t slice_end = clock();
  slicing_time = double(slice_end - slice_begin)/CLOCKS_PER_SEC;
  /*End-Slicing*/

  /*Changing structure for CGAL library: */
  vector<Segment_2> S[k];
  for (int p = 0; p < k; p++) {
    for (auto it = segs[p].begin(), itEnd = segs[p].end(); it != itEnd; ++it) {
      LineSegment s = *it; /*Current segment.*/
      S[p].push_back(Segment_2(Point_2(s.v[0].x, s.v[0].y), Point_2(s.v[1].x, s.v[1].y)));
    }
  } 

  if (chaining) {
    /*Loop-Closure:*/
    clock_t contour_begin = clock();
    for (size_t p = 0; p < k; p++) {
      if (!S[p].empty()) {
         ParkLoopClosure (S[p], polygons);
      }
    }
    clock_t contour_end = clock();
    loopclosure_time = double(contour_end - contour_begin)/CLOCKS_PER_SEC;
    /*End-Loop-Closure*/
  }
  else {
    export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
  }
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*******************************  INPUT/OUTPUT and Auxiliar routine **********************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

/*-----------------------------------------------------------------------*/
/*Rounds x to an integer multiple of {eps}. If {mod} > 1,
  the factor will be congruent to {rem} modulo {mod}. */
float xround (float x, double eps, int mod, int rem) {
  double y = round((double)x/(mod * eps));
  double z = (y * mod + rem) * eps;
  return (float)z;
}

/*-----------------------------------------------------------------------*/
/*Rounds {x,y,z} to an even multiple of {eps}. */
v3 v3_round (float x, float y, float z, double eps) {
   v3 p;
   p.x = xround(x, eps, 2, 0);
   p.y = xround(y, eps, 2, 0);
   p.z = xround(z, eps, 2, 0);
   return p;
}

/*-----------------------------------------------------------------------*/
/*Returns true {t} has two or more coicident vertices.*/
bool degenerate (Triangle t) {
   if (t.v[0].distTo(t.v[1]) < 0.000001) { return true; }
   if (t.v[1].distTo(t.v[2]) < 0.000001) { return true; }
   if (t.v[2].distTo(t.v[0]) < 0.000001) { return true; }
   return false;
}

/*-----------------------------------------------------------------------*/
Triangle make_triangle ( 
   float n0, float n1, float n2,
   float f0, float f1, float f2,
   float f3, float f4, float f5,
   float f6, float f7, float f8,
   const char *rotate, double eps
 ) {

    if (strcmp(rotate,"true") == 0) {
      return Triangle (v3(n0, n1, n2), v3_round(f0, f2, f1, eps), v3_round(f3, f5, f4, eps), v3_round(f6, f8, f7, eps));
    }
    else {
      return Triangle (v3(n0, n1, n2), v3_round(f0, f1, f2, eps), v3_round(f3, f4, f5, eps), v3_round(f6, f7, f8, eps));
    }
}

/*-----------------------------------------------------------------------*/
/* Read the given STL file name (ascii or binary is set using ‘isBinaryFormat’)
   and generate a Triangle Mesh object in output parameter ‘mesh’. */
int stlToMeshInMemory (const char *stlFile, TriangleMesh *mesh, bool isBinaryFormat, const char *rotate, double eps) {

  int ndegenerated = 0;

  if (!isBinaryFormat) {

    ifstream in(stlFile);

    if (!in.good()) {
      return 1;
    }

    std::string s0, s1;

    float n0, n1, n2;
    float f0, f1, f2;
    float f3, f4, f5;
    float f6, f7, f8;

    while (!in.eof()) {
      in >> s0; /*s0 can be facet or solid!*/
      if (s0 == "facet") {
        in >> s1 >> n0 >> n1 >> n2; /* normal x y z. */
        in >> s0 >> s1;             /* loop. */
        in >> s0 >> f0 >> f1 >> f2; /* vertex x y z. */
        in >> s0 >> f3 >> f4 >> f5; /* vertex x y z. */
        in >> s0 >> f6 >> f7 >> f8; /* vertex x y z. */
        in >> s0;                   /* endloop. */
        in >> s0;                   /* endfacet.*/
        Triangle t = make_triangle (n0, n1, n2, f0, f2, f1, f3, f5, f4, f6, f8, f7, rotate, eps);
        if (!degenerate(t)) {
           mesh->push_back (t);
        }
        else {
          ndegenerated++;
        }
      }
      else if (s0 == "endsolid") {
         break;
      }
    }
    in.close();
  }
  else {
    FILE *f = fopen (stlFile, "rb");
    if (!f) {
      return 1;
    }
    char title[80];
    int nFaces;
    int err;
    err = fread (title, 80, 1, f);
    err = fread ((void*)&nFaces, 4, 1, f);
    float v[12]; /* normal = 3, vertices = 9 (12) */
    unsigned short uint16;
    /* Every Face is 50 Bytes: Normal(3*float), Vertices(9*float), 2 Bytes Spacer */
    for (size_t i=0; i<nFaces; ++i) {
      for (size_t j=0; j<12; ++j) {
        err = fread((void*)&v[j], sizeof(float), 1, f);
      }
      err = fread((void*)&uint16, sizeof(unsigned short), 1, f); // spacer between successive faces
      Triangle t = make_triangle (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], rotate, eps);
      if (!degenerate(t)) {
         mesh->push_back (t);
      }
      else {
        ndegenerated++;
      }
    }
    fclose(f);
  }
  cout << "number of degenerated triangles = " << ndegenerated << endl;
  return 0;
}

/*-----------------------------------------------------------------------*/
/*Compute uniform and adaptive z-plane coordinates!*/
vector<float> compute_planes (const TriangleMesh *mesh, float max_thickness, char *adaptive, double eps, float *delta) {

  bool rounding = true; /*To avoid that the Z-coordinates of all planes are distinct from the Z-coordinates of all vertices.*/

  /* Vector to keep the plane coordinates: */
  vector<float> Planes;

  /* Assuming the model as a 3D axis-aligned bounding-box: */
  double model_zmax = std::max(mesh->getUpperRightVertex().z, mesh->meshAABBSize().z);

  double model_zmin = mesh->getBottomLeftVertex().z;

  if (strcmp(adaptive, "false") == 0)  { /*Uniform slicing: */

    double spacing = (rounding ? xround (max_thickness, eps, 2, 0) : max_thickness); /*Plane spacing even multiple of {eps}*/

    double P0 = xround (model_zmin - spacing, eps, 2, 1); /*First plane odd multiple of {eps}.*/

    int no_planes = 1 + (int)((model_zmax + spacing - P0)/spacing); /* Number of planes: */

    cout << "eps = " << eps << endl;
    cout << "max thickness = " << max_thickness << endl;
    cout << "rounded plane spacing spacing = " << spacing << endl;
    cout << "model zmin = " << model_zmin << ", model zmax = " << model_zmax << ", first plane Z = " << P0 << ", number of planes = " << no_planes << endl;

    for (size_t i = 0; i < no_planes; i++) {
      /* Building the vector with the slice z coordinates: */
      float Pi = (float)(P0 + i * spacing);
      if ((Pi > model_zmin) && (Pi < model_zmax)) {
          Planes.push_back ((float)(P0 + i * spacing));
      }
    }
    *delta = (float)(spacing);
  }
  else { /*Adaptive slicing z-planes: */

    float zplane = 0.0;
    float min_thickness = 0.016;
    Planes.push_back (model_zmin + zplane);

    while ((model_zmin + zplane) <= model_zmax) {
      double vrandom = min_thickness + (max_thickness - min_thickness) * (rand() / (double)RAND_MAX);
      double coordinate = xround (model_zmin + zplane + vrandom, eps, 2, 1);
      if (coordinate >= model_zmax) { break; }
      Planes.push_back (coordinate);
      zplane += vrandom;
    }
  }
  return Planes;
}
/*-----------------------------------------------------------------------*/
void fill_colors (rgb colors[], int nrgb) {
  assert (nrgb == 8);
  colors[0] = {  0,   0,   0}; /*Black*/
  colors[1] = {255,   0,   0}; /*Red*/
  colors[2] = {128,   0, 128}; /*Purple*/
  colors[3] = {  0, 128,   0}; /*Green*/
  colors[4] = {  0,   0, 255}; /*Blue*/
  colors[5] = {  0, 255, 255}; /*Cyan*/
  colors[6] = {128, 128,   0}; /*Olive*/
  colors[7] = {128, 128, 128}; /*Gray*/
}

/*-----------------------------------------------------------------------*/
/*SVG header information:*/
void add_svg_information (FILE *file) {
  fprintf (file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf (file, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf (file, "<svg preserveAspectRatio=\"xMidYMid meet\" width=\"1024\" height=\"768\" viewBox=\"0 0 1024 768\""); 
  fprintf (file, " xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"");
  fprintf (file, " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\""); 
  fprintf (file, " xmlns:cc=\"http://web.resource.org/cc/\">\n");
  fprintf (file, " <defs>\n");
  fprintf (file, "  <marker id=\"arrow\" markerWidth=\"10\" markerHeight=\"10\" refx=\"0\" refy=\"3\" orient=\"auto\" markerUnits=\"strokeWidth\" viewBox=\"0 0 20 20\">\n");
  fprintf (file, "    <path d=\"M0,0 L0,6 L9,3 z\" fill=\"#f00\" />\n");
  fprintf (file, "  </marker>\n");
  fprintf (file, "</defs>\n");
}


/*-----------------------------------------------------------------------*/
/*Output a SVG file with the whole 3D model*/
void export_svg (char filename[], const vector<vector<v3> > &Polygons, int nslices, int nsegments, bool video, vector<bool> orientation) {
  glm::vec3 fromEuler (0.0f, 0.0f, 60.0f);
  glm::quat quaternion (DEG_TO_RAD(fromEuler));
  glm::vec3 toEuler = glm::eulerAngles(quaternion);
  float angle = glm::angle(quaternion);
  glm::vec3 axis = glm::axis(quaternion);
  glm::mat4 View = glm::rotate(glm::mat4(1.0), angle, axis);

  //float zoom = 0.6f; /*01.liver*/
  //float zoom = 0.2f; /*02.femur*/
  float zoom = 0.25f;  /*04.demon*/
  glm::mat4 Projection = glm::perspective (zoom, 4.0f / 3.0f, 0.1f, 100.f);
  int shift_x = +400;
  int shift_y = +400;

  glm::mat4 Model = glm::lookAt (
        glm::vec3(1, 1, 1),    /* Eye point (where am I?) */
        glm::vec3(0, 0, 0),    /* Look at Center point (What am I looking at?) */
        glm::vec3(0, 1, 0));   /* UP vector of the camera (Where is my UP vector?) */
  glm::mat4 MVP = Projection * View * Model;
  FILE *file = fopen (filename, "w");
  add_svg_information (file);

  /*Writing previous slices:*/
  int i;
  for (i = 0; i < nslices; i++) {
    for (size_t index = 0; index < Polygons[i].size(); index++) {
      v3 p0 = Polygons[i].at(index);
      if (index < (Polygons[i].size() - 1)) {
        v3 p1 = Polygons[i].at(index + 1);
        p0.transform(MVP);
        p1.transform(MVP);
        if (orientation.at(i)) {
          fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(255,0,0)\"/>\n",
                   p0.x + shift_x, p0.y + shift_y, p1.x + shift_x,p1.y + shift_y);
        }
        else {
          fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,0)\"/>\n",
                   p0.x + shift_x, p0.y + shift_y, p1.x + shift_x,p1.y + shift_y);
        }
      }
    }
  }
  if (video) {
    /*Writing current slice:*/
    for (int j = i; j <= nslices; j++) {
      for (int index = 0; index < nsegments; index++) {
        v3 p0 = Polygons[j].at(index);
        if (index < (Polygons[j].size() - 1)) {
          v3 p1 = Polygons[j].at(index + 1);
          p0.transform(MVP);
          p1.transform(MVP);
          if (orientation.at(i)) {
             fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(255,0,0)\"/>\n",
                      p0.x + shift_x, p0.y + shift_y, p1.x + shift_x,p1.y + shift_y);
          }
          else {
             fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(0,0,0)\"/>\n",
                      p0.x + shift_x, p0.y + shift_y, p1.x + shift_x,p1.y + shift_y);
          }
        }
      }
    }
  }
  fprintf (file,"</svg>\n");
  fclose(file);
}

/*-----------------------------------------------------------------------*/
/*Output a SVG file with each model layer (list of oriented points in {Polygons})*/
void export_svg_2d (string fileName, vector<contour> polygons[], int nplanes, const v3 &aabbSize) {

  int nrgb = 8;
  rgb colors[nrgb];
  fill_colors (colors, nrgb);

  FILE *file = NULL;

  float dx = 0.0, dy = 0.0;

  file = fopen (fileName.c_str(), "w");
  printf("\n\nwriting output file: %s\n", fileName.c_str());

  if (!file) { exit(1); }

  add_svg_information (file);
    
  for (int p = 0; p < nplanes; p++) { 
    vector<contour> P = polygons[p];
    const size_t k = P.size();
    const size_t slicePerRow = (size_t)sqrt((float)nplanes)*2;
    for (size_t i = 0; i < k; ++i) {
      for (size_t index = 0; index < P[i].P.size(); index++) {
        const v3 &p0 = P[i].P.at(index);
        dx = (float)(p % slicePerRow)*(aabbSize.x * 1.05f);
        dy = (float)(p / slicePerRow)*(aabbSize.y * 1.05f);
        if (index < (P[i].P.size() - 1)) {
          const v3 &p1 = P[i].P.at(index + 1);
          //if (P[i].external) {
          //if (true) {
          if (P[i].clockwise) {
            //fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\" marker-end=\"url(#arrow)\"/> \n",
            fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\"/>\n",
            dx + p0.x, dy + p0.y, dx + p1.x, dy + p1.y, colors[0].r, colors[0].g, colors[0].b);
            //dx + p0.x, dy + p0.y, dx + p1.x, dy + p1.y, colors[i % nrgb].r, colors[i % nrgb].g, colors[i % nrgb].b);
          }
          else {
            //fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\" marker-end=\"url(#arrow)\"/>\n",
            fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\"/>\n",
            dx + p0.x, dy + p0.y, dx + p1.x, dy + p1.y, colors[1].r, colors[1].g, colors[1].b); 
            //dx + p0.x, dy + p0.y, dx + p1.x, dy + p1.y, colors[i % nrgb].r, colors[i % nrgb].g, colors[i % nrgb].b);
          }
        }
      }
    }
  }
  fprintf (file,"</svg>\n");
  fclose (file);
}

/*Output a SVG file without chaining (list of unoriented line segments).*/
void export_svg_no_chaining (string fileName, vector<LineSegment> Segments[], int k, const v3 &aabbSize) {

  int nrgb = 8;
  rgb colors[nrgb];
  fill_colors (colors, nrgb);

  FILE *file = NULL;

  float dx = 0.0, dy = 0.0;

  file = fopen (fileName.c_str(), "w");
  printf("\n\nwriting output file: %s\n", fileName.c_str());

  if (!file) { exit(1); }

  add_svg_information (file);

  const size_t slicePerRow = (size_t)sqrt((float)k)*2;

  for (size_t i = 0; i < k; ++i) {
    dx = (float)(i%slicePerRow)*(aabbSize.x * 1.05f);
    dy = (float)(i/slicePerRow)*(aabbSize.y * 1.05f);
    for (const LineSegment &ls : Segments[i]) {
    fprintf (file, "   <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"rgb(%d,%d,%d)\"/>\n",
             dx + ls.v[0].x, dy + ls.v[0].y, dx + ls.v[1].x, dy + ls.v[1].y,
             colors[i % nrgb].r, colors[i % nrgb].g, colors[i % nrgb].b);
    }
    Segments[i].clear(); 
  }
  fprintf (file,"</svg>\n");
  fclose (file);
}

/*-----------------------------------------------------------------------*/
void export_video (vector<contour> polygons[], int nplanes) {

  vector<vector<v3>> P;
  vector<bool> orientation;
  for (int k = 0; k < nplanes; k++) {
    const size_t ncontours = polygons[k].size();
    for (int c = 0; c < ncontours; c++) {
      vector<v3> Pc = polygons[k].at(c).P;
      orientation.push_back(polygons[k].at(c).clockwise);
      P.push_back(Pc);
    }
  }

  for (int k = 0; k < P.size(); k++) {
    const size_t nsegments = P.at(k).size();
    for (int j = 0; j < nsegments; j++) {
      char filename[256];
      sprintf(filename, "./video/%05d_%05d.svg", k, j);
      printf("writing output file: %s\n", filename);
      export_svg (filename, P, k, j, true, orientation);
      printf("... done\n");
    }
  }
  printf("... done\n");
}

/*-----------------------------------------------------------------------*/
void export_svg_3d (vector<contour> polygons[], int nplanes) {

  vector<vector<v3>> P;
  vector<bool> orientation;
  for (int k = 0; k < nplanes; k++) {
    const size_t ncontours = polygons[k].size();
    for (int c = 0; c < ncontours; c++) {
      vector<v3> Pc = polygons[k].at(c).P;
      orientation.push_back(polygons[k].at(c).clockwise);
      P.push_back(Pc);
    }
  }

  char filename[256];
  sprintf(filename, "out_3d.svg");
  printf("\n\nwriting output file: %s\n", filename);
  export_svg (filename, P, P.size(), 0, false, orientation);
  printf("... done\n\n\n");
}

/*-----------------------------------------------------------------------*/
int checkASCIIFile (const char *fileName) {
  string line1, line2;
  ifstream input(fileName);
  if (!input.good()) {
    return -1;
  }
  getline(input, line1);
  getline(input, line2);
  if (line1.find("solid")!=string::npos && line2.find("facet")!=string::npos) {
    return FILE_STL_ASCII;
  }
  if (line1.find("xml")!=string::npos && line2.find("amf")!=string::npos) {
    return FILE_AMF_ASCII;
  }
  return FILE_STL_BIN;
}

/*----------------------------Main procedure---------------------------------*/
int main (int argc, char **argv) {

  bool chaining = true;

  double eps = 0.004;

  /*Total time:*/
  clock_t begin = clock();

  char *model;

  if (strcmp(argv[2], "-model") == 0) {
     model = argv[3];
  }

  float max_thickness, delta;

  if (strcmp(argv[4], "-thickness") == 0) {
    max_thickness = atof(argv[5]);
  }  
  else {
    printf("Error: specify the slicing spacing in mm (thickness)!!!\n");
  }

  char *adaptive;

  if (strcmp(argv[6], "-adaptive") == 0) {
    adaptive = argv[7];
  }  

  char *write_option = argv[8];

  char *rotate;

  if (strcmp(argv[9], "-rotate") == 0) {
     rotate = argv[10];
  }

  bool orienting;

  if (strcmp(argv[11], "-orienting_contours") == 0) {
    if (strcmp(argv[12],"true") == 0) {
      orienting = true;
    }
    else {
      orienting = false;
    }
  }

  TriangleMesh mesh;
    
  switch (checkASCIIFile(model)) {
    case FILE_STL_ASCII:
      if (stlToMeshInMemory (model, &mesh, false, rotate, eps) != 0)
        return 1;
      break;
    case FILE_STL_BIN:
      if (stlToMeshInMemory (model, &mesh, true, rotate, eps) != 0)
        return 1;
      break;
    default:
      cerr << "Unexpected error" << endl;
      return 1;
  }

  std::string path = model;
    
  std::string lastFileName;
    
  size_t pos = path.find_last_of("/");
    
  if (pos != std::string::npos)
    lastFileName.assign(path.begin() + pos + 1, path.end());
  else
    lastFileName = path;
    
  vector<float> P = compute_planes (&mesh, max_thickness, adaptive, eps, &delta);
  
  int nplanes = P.size();

  vector<contour> polygons[nplanes];
 
  bool srt = false;
  
  if (strcmp(write_option, "-No") == 0) {
     chaining = false;
  }

  if (strcmp(argv[1], "-Trivial") == 0) {
    TrivialSlicing (&mesh, P, polygons, chaining, orienting);
  } 
  else if (strcmp(argv[1], "-Park") == 0) {
    ParkSlicing (&mesh, P, polygons, chaining);
  }
  else if (strcmp(argv[1], "-Incremental") == 0) {
    if (strcmp(adaptive, "true") == 0) {
      delta = 0.0;
    }
    IncrementalSlicing (&mesh, P, delta, srt, polygons, chaining, orienting);
  }
  
  clock_t end = clock();

  double total_time = double(end - begin)/CLOCKS_PER_SEC;

  cout << argv[1] << ", " << lastFileName << ", thickness = " << delta; 
  cout << ", #T = " << mesh.size() << ", #P = " << nplanes; cout << ", #S = " << intersections;
  cout << ", " << slicing_time << "s (slicing), " << loopclosure_time << "s (polygon assembling), " << total_time << "s (total)";
  cout << ", *K: " << ((double)intersections/(double)mesh.size()) << endl;

  if (chaining) {
    if (strcmp(write_option, "-2D") == 0) {
      export_svg_2d ("out_2d.svg", polygons, nplanes, mesh.meshAABBSize());
    }
    else if (strcmp(write_option, "-3D") == 0) {
      export_svg_3d (polygons, nplanes);
    }  
    else if (strcmp(write_option, "-video") == 0) {
      export_video (polygons, nplanes);
    }
  }
  
  return 0;
}

