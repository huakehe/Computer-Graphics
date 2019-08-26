//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_BOUNDS3_H
#define RAYTRACING_BOUNDS3_H
#include <limits>
#include "Vector.hpp"
#include "Ray.hpp"
#include <algorithm>

class Bounds3{
public:
    Vector3f pMin,pMax; // two points to specify the bounding box
    Bounds3(){
        double minNum = std::numeric_limits<double>::lowest();
        double maxNum = std::numeric_limits<double>::max();
        pMax = Vector3f(minNum,minNum,minNum);
        pMin = Vector3f(maxNum,maxNum,maxNum);
    }
    Bounds3(const Vector3f p):pMin(p), pMax(p){}
    Bounds3(const Vector3f p1, const Vector3f p2){
        pMin = Vector3f(fmin(p1.x,p2.x),fmin(p1.y,p2.y),fmin(p1.z,p2.z));
        pMax = Vector3f(fmax(p1.x,p2.x),fmax(p1.y,p2.y),fmax(p1.z,p2.z));
    }

    Vector3f Diagonal()const {return pMax-pMin;}
    int maxExtent()const{
        Vector3f d = Diagonal();
        if(d.x>d.y && d.x>d.z)
            return 0;
        else if (d.y>d.z)
            return 1;
        else
            return 2;
    }

    double SurfaceArea()const{Vector3f d = Diagonal();return 2*(d.x*d.y+d.x*d.z+d.y*d.z);}

    Vector3f Centroid(){return 0.5*pMin + 0.5*pMax;}
    Bounds3 Intersect(const Bounds3 &b){
        return Bounds3(Vector3f(fmax(pMin.x,b.pMin.x),
                               fmax(pMin.y,b.pMin.y),
                               fmax(pMin.z,b.pMin.z)),
                       Vector3f(fmin(pMax.x,b.pMax.x),
                               fmin(pMax.y,b.pMax.y),
                               fmin(pMax.z,b.pMax.z))
        );
    }

    Vector3f Offset(const Vector3f &p) const {
        Vector3f o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
        return o;
    }

    bool Overlaps(const Bounds3 &b1, const Bounds3 &b2) {
        bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
        bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
        bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
        return (x && y && z);
    }

    bool Inside(const Vector3f &p, const Bounds3 &b) {
        return (p.x >= b.pMin.x && p.x <= b.pMax.x &&
                p.y >= b.pMin.y && p.y <= b.pMax.y &&
                p.z >= b.pMin.z && p.z <= b.pMax.z);
    }
    inline const Vector3f& operator[](int i) const {
        return (i == 0) ? pMin : pMax;
    }

    inline bool IntersectP(const Ray &ray, float* hitt0, float* hitt1)const;
    inline bool IntersectP(const Ray &ray, const Vector3f &invDir, const int dirisNeg[3])const;
};



inline bool Bounds3::IntersectP(const Ray &ray, float* hitt0, float* hitt1)const{
    float t0 = 0, t1 = ray.t_max;
    for(int i=0;i<3;i++){
        float invRay = ray.direction_inv[i];
        float tNear = (pMin[i]-ray.origin[i])*invRay;
        float tFar = (pMax[i]-ray.origin[i])*invRay;
        if(tNear>tFar)std::swap(tNear,tFar);
        t0 = tNear > t0 ? tNear:t0;
        t1 = tFar  < t1 ? tFar :t1;
        if(t0>t1)return false;
    }

    if(hitt0) *hitt0 = t0;
    if(hitt1) *hitt1 = t1;
    return true;
}

inline bool Bounds3::IntersectP(const Ray &ray, const Vector3f &invDir, const int dirIsNeg[3])const{
    //TODO test if ray bound intersects
    float t0 = 0, t1 = ray.t_max;
    for(int i=0;i<3;i++){
        float invRay = ray.direction_inv[i];
        float tNear = (pMin[i]-ray.origin[i])*invRay;
        float tFar = (pMax[i]-ray.origin[i])*invRay;
        if(tNear>tFar)std::swap(tNear,tFar);
        t0 = tNear > t0 ? tNear:t0;
        t1 = tFar  < t1 ? tFar :t1;
        if(t0>t1)return false;
    }

    return true;
    /*
    double txmin = (pMin.x-ray.origin.x)*invDir.x;
    double txmax = (pMax.x-ray.origin.x)*invDir.x;
    double tymin = (pMin.y-ray.origin.y)*invDir.y;
    double tymax = (pMax.y-ray.origin.y)*invDir.y;
    double tzmin = (pMin.z-ray.origin.z)*invDir.z;
    double tzmax = (pMax.z-ray.origin.z)*invDir.z;
    
    
    double tenter,texit;
    tenter = std::max(txmin,tymin);
    tenter = std::max(tenter,tzmin);
    texit = std::min(txmax,tymax);
    texit = std::min(texit,tzmax);
    if(texit<0) return false;
    if(tenter<0) return true;
    if(tenter>=texit) return false;
    else{
        return true;
    }*/
    
    /*
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    
    if (invDir.x >= 0) {
        tmin = (pMin.x - ray.origin.x) * invDir.x;
        tmax = (pMax.x - ray.origin.x) * invDir.x;
    }
    else {
        tmin = (pMax.x - ray.origin.x) * invDir.x;
        tmax = (pMin.x - ray.origin.x) * invDir.x;
    }
    
    if (invDir.y >= 0) {
        tmin = (pMin.y - ray.origin.y) * invDir.y;
        tmax = (pMax.y - ray.origin.y) * invDir.y;
    }
    else {
        tmin = (pMax.y - ray.origin.y) * invDir.y;
        tmax = (pMin.y - ray.origin.y) * invDir.y;
    }
    
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    
    if (invDir.z >= 0) {
        tmin = (pMin.z - ray.origin.z) * invDir.z;
        tmax = (pMax.z - ray.origin.z) * invDir.z;
    }
    else {
        tmin = (pMax.z - ray.origin.z) * invDir.z;
        tmax = (pMin.z - ray.origin.z) * invDir.z;
    }
    
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    
    return true;*/
    
     /*
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    Vector3f bounds[2];
    bounds[0] = pMin;
    bounds[1] = pMax;
    
    tmin = (bounds[dirIsNeg[0]].x - ray.origin.x) * invDir.x;
    tmax = (bounds[1-dirIsNeg[0]].x - ray.origin.x) * invDir.x;
    
    tymin = (bounds[dirIsNeg[1]].y - ray.origin.y) * invDir.y;
    tymax = (bounds[1-dirIsNeg[1]].y - ray.origin.y) * invDir.y;
    
    
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    
    tmin = (bounds[dirIsNeg[2]].z - ray.origin.z) * invDir.z;
    tmax = (bounds[1-dirIsNeg[2]].z - ray.origin.z) * invDir.z;
    
    
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    
    return true;*/
}



inline Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2) {
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b1.pMin, b2.pMin);
    ret.pMax = Vector3f::Max(b1.pMax, b2.pMax);
    return ret;
}

inline Bounds3 Union(const Bounds3 &b,const Vector3f &p){
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b.pMin, p);
    ret.pMax = Vector3f::Max(b.pMax, p);
    return ret;
}

#endif //RAYTRACING_BOUNDS3_H
