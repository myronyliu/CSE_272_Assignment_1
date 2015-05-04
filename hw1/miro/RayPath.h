#ifndef CSE168_RAYPATH_H_INCLUDED
#define CSE168_RAYPATH_H_INCLUDED

#include "Ray.h"
#include <vector>
using namespace std;

class RayPath
{
public:
    Ray rayInit;
    vector<Ray> rays;
    vector<HitInfo> hits;
    vector<float> probs;
};

#endif // CSE168_RAYPATH_H_INCLUDED