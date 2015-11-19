#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform
#include <OpenGL/gl.h>
#include "float2.h"
#include <float.h>
#include <OpenGL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GLUT/glut.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>

// simple material class, with object color, and headlight shading
class Material
{
public:
    virtual float3 getColor(
                            float3 position,
                            float3 normal,
                            float3 viewDir)
    {
        return normal;
    }
    virtual float3 shade(
                         float3 normal,
                         float3 viewDir,
                         float3 lightDir,
                         float3 lightPowerDensity,
                         float3 position)
    {
        return normal;
    };
};

class LightSource
{
public:
    virtual float3 getPowerDensityAt ( float3 x )=0;
    virtual float3 getLightDirAt     ( float3 x )=0;
    virtual float  getDistanceFrom   ( float3 x )=0;
};

// class DirectionalLight : public LightSource
class DirectionalLight: public LightSource
{
    float3 direction;
    float3 powerDensity;
public:
    DirectionalLight(float3 direction, float3 powerDensity):
    direction(direction),
    powerDensity(powerDensity)
    {
    }
    float3 getPowerDensityAt(float3 x) {
        return powerDensity;
    }
    
    float3 getLightDirAt(float3 x) {
        return direction;
    }
    float getDistanceFrom(float3 x) {
        return 10000;
    }
    
};
// class PointLight : public LightSource
class PointLight: public LightSource{
    float3 position;
    float3 powerDensity;
public:
    PointLight(float3 position, float3 powerDensity):
    position(position),
    powerDensity(powerDensity)
    {
    }
    
    float3 getPowerDensityAt(float3 x) {
        float distance = position.distance(x);
        float3 newPowerDensity = powerDensity * (1/(distance*distance));
        return newPowerDensity;
    }
    
    float3 getLightDirAt(float3 x) {
        float3 direction = position - x;
        direction.normalize();
        return direction;
    }
    
    float getDistanceFrom(float3 x) {
        float distance = position.distance(x);
        return distance;
    }
};


// Skeletal camera class.
class Camera
{
    float3 eye;		//< world space camera position
    float3 lookAt;	//< center of window in world space
    float3 right;	//< vector from window center to window right-mid (in world space)
    float3 up;		//< vector from window center to window top-mid (in world space)
    
public:
    Camera()
    {
        eye = float3(0, 0, 3);
        lookAt = float3(0, 0, 2);
        right = float3(1, 0, 0);
        up = float3(0, 1, 0);
    }
    float3 getEye()
    {
        return eye;
    }
    // compute ray through pixel at normalized device coordinates
    float3 rayDirFromNdc(const float2 ndc) {
        return (lookAt - eye
                + right * ndc.x
                + up    * ndc.y
                ).normalize();
    }
};

class HeadlightMaterial : public Material {
    float3 frontFaceColor;
    float3 backFaceColor;
public:
    HeadlightMaterial(float3 frontfaceColor,
                      float3 backfaceColor  ):
    frontFaceColor(frontFaceColor),
    backFaceColor(backFaceColor){}
    HeadlightMaterial():
    frontFaceColor(float3::random()),
    backFaceColor(float3::random())
    {}
    float3 getColor(
                            float3 position,
                            float3 normal,
                            float3 viewDir) {
        //implement headlight shading formula here
        float3 newColor;
        if (viewDir.dot(normal) < 0) {
            newColor.x = -viewDir.dot(normal)*backFaceColor.x;
            newColor.y = -viewDir.dot(normal)*backFaceColor.y;
            newColor.z = -viewDir.dot(normal)*backFaceColor.z;
        }
        else {
            newColor.x = viewDir.dot(normal)*frontFaceColor.x;
            newColor.y = viewDir.dot(normal)*frontFaceColor.y;
            newColor.z = viewDir.dot(normal)*frontFaceColor.z;
        }
        
        return newColor;
        
    }
    //TODO: fix this!!!!
    float3 shade(
                 float3 normal,
                 float3 viewDir,
                 float3 lightDir,
                 float3 lightPowerDensity,
                 float3 position){
        return frontFaceColor;
    };
};


class Metal : public Material {
    float3 r0;
public:
    Metal(float3  refractiveIndex, float3  extinctionCoefficient){
        float3 rim1 = refractiveIndex - float3(1,1,1);
        float3 rip1 = refractiveIndex + float3(1,1,1);
        float3 k2 = extinctionCoefficient * extinctionCoefficient;
        r0 = (rim1*rim1 + k2) / (rip1*rip1 + k2);
    }
    
    struct Event{
        float3 reflectionDir;
        float3 reflectance;
    };
    
    Event evaluateEvent(float3 inDir, float3 normal) {
        Event e;
        float cosa = -normal.dot(inDir);
        float3 perp = -normal * cosa;
        float3 parallel = inDir - perp;
        e.reflectionDir = parallel - perp;
        
        e.reflectance = r0 + (float3(1,1,1)-r0) * pow(1 - cosa, 5);
        
        return e; }
};

class ProceduralTextureMetal : public Material {
    float3 r0;
public:
    ProceduralTextureMetal(float3  refractiveIndex, float3  extinctionCoefficient){
        float3 rim1 = refractiveIndex - float3(1,1,1);
        float3 rip1 = refractiveIndex + float3(1,1,1);
        float3 k2 = extinctionCoefficient * extinctionCoefficient;
        r0 = (rim1*rim1 + k2) / (rip1*rip1 + k2);
    }
    
    struct Event{
        float3 reflectionDir;
        float3 reflectance;
    };
    
    //gradient of noise function
    float3 snoiseGrad(float3 r)
    {
        unsigned int x = 0x0625DF73;
        unsigned int y = 0xD1B84B45;
        unsigned int z = 0x152AD8D0;
        float3 f = float3(0, 0, 0);
        for(int i=0; i<32; i++)
        {
            float3 s( x/(float)0xffffffff,
                     y/(float)0xffffffff,
                     z/(float)0xffffffff);
            f += s * cos(s.dot(r));
            x = x << 1 | x >> 31;
            y = y << 1 | y >> 31;
            z = z << 1 | z >> 31;
        }
        return f * (1.0 / 64.0);
    }
    
    //modify this to take in position
    Event evaluateEvent(float3 inDir, float3 normal, float3 position) {
        float3 sNoiseGrad = snoiseGrad(position).normalize();
        float freq = 10;
        float amplitude = 1;
        float3 perturbedNormal = normal + snoiseGrad(position * freq) * amplitude;
        perturbedNormal.normalize();
        
        Event e;
        float cosa = -perturbedNormal.dot(inDir);
        float3 perp = -perturbedNormal * cosa;
        float3 parallel = inDir - perp;
        e.reflectionDir = parallel - perp;
        
        e.reflectance = r0 + (float3(1,1,1)-r0) * pow(1 - cosa, 5);
        
        return e; }
    
};

class PhongBlinn : public Material {
    float3 ks;
    float3 kd;
    float shininess;
public:
    PhongBlinn(float3 ks, float3 kd, float shininess):
    ks(ks),
    kd(kd),
    shininess(shininess) {}
    float3 shade( float3 normal, float3 viewDir,
                 float3 lightDir, float3 lightPowerDensity, float3 position)
    {
        float3 halfway =
        (viewDir + lightDir).normalize();
        float cosDelta = normal.dot(halfway);
        float cosTheta = normal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        
        float3 variable =  kd * lightPowerDensity * cosTheta;
        
        // return variable;
        
        if(cosDelta < 0) return variable;
        return lightPowerDensity * ks
        * pow(cosDelta, shininess) + variable;
        
    }
};

class ProceduralTexturePhongBlinn : public Material{
    float3 ks;
    float3 kd;
    float shininess;
public:
    ProceduralTexturePhongBlinn(float3 ks, float3 kd, float shininess):
    ks(ks),
    kd(kd),
    shininess(shininess)
    {}
    
    //noise function
    float snoise(float3 r) {
        unsigned int x = 0x0625DF73;
        unsigned int y = 0xD1B84B45;
        unsigned int z = 0x152AD8D0;
        float f = 0;
        for(int i=0; i<32; i++) {
            float3 s(	x/(float)0xffffffff,
                     y/(float)0xffffffff,
                     z/(float)0xffffffff);
            f += sin(s.dot(r));
            x = x << 1 | x >> 31;
            y = y << 1 | y >> 31;
            z = z << 1 | z >> 31;
        }
        return f / 64.0 + 0.5;
    }
    
    //gradient of noise function
    float3 snoiseGrad(float3 r)
    {
        unsigned int x = 0x0625DF73;
        unsigned int y = 0xD1B84B45;
        unsigned int z = 0x152AD8D0;
        float3 f = float3(0, 0, 0);
        for(int i=0; i<32; i++)
        {
            float3 s( x/(float)0xffffffff,
                     y/(float)0xffffffff,
                     z/(float)0xffffffff);
            f += s * cos(s.dot(r));
            x = x << 1 | x >> 31;
            y = y << 1 | y >> 31;
            z = z << 1 | z >> 31;
        }
        return f * (1.0 / 64.0);
    }
    
    float3 shade(
                 float3 normal,
                 float3 viewDir,
                 float3 lightDir,
                 float3 lightPowerDensity,
                 float3 position)
    {
        //return based on position
        //normal mapping: use sNoisegrad function to return gradient
        float3 sNoiseGrad = snoiseGrad(position).normalize();
        //add this color to the normal
        float freq = 110;
        float amplitude = 5;
        float3 perturbedNormal = normal + snoiseGrad(position * freq) * amplitude;
        perturbedNormal.normalize();
        float3 halfway =
        (viewDir + lightDir).normalize();
        float cosDelta = perturbedNormal.dot(halfway);
        float cosTheta = perturbedNormal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        
        float3 variable =  kd * lightPowerDensity * cosTheta;
        
        // return variable;
        
        if(cosDelta < 0) return variable;
        return (lightPowerDensity * ks
        * pow(cosDelta, shininess) + variable);

        
    }
    
};


class DiffusePlane : public Material
{
    float3 kd;
    //for chessboard: kd depends on where hit position is
    
public:
    
    float3 getColor(
                    float3 position,
                    float3 normal,
                    float3 viewDir)
    {
        return normal;
    }
    
    float3 shade(
                 float3 normal,
                 float3 viewDir,
                 float3 lightDir,
                 float3 lightPowerDensity,
                 float3 position)
    {
        //add some large positive integer before converting
        
        if (((int)(position.x + 100)+ (int)(position.y + 100) + (int)(position.z + 100) ) % 2 ==0) {
            kd = float3(0.0, 0.0, 0.0);
        }
        else {
            kd = float3(1.0, 1.0, 1.0);
        }
        //        return kd * lightPowerDensity * cosTheta;
        return kd;
    }
};


// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
    Hit()
    {
        t = -1;
    }
    float t;				//< Ray paramter at intersection. Negative means no valid intersection.
    float3 position;		//< Intersection coordinates.
    float3 normal;			//< Surface normal at intersection.
    Material* material;		//< Material of intersected surface.
};

// Object abstract base class.
class Intersectable
{
protected:
    Material* material;
public:
    Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
    
    
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.
class QuadraticRoots
{
public:
    float t1;
    float t2;
    // Solves the quadratic a*a*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and set members t1 and t2 to store the roots.
    QuadraticRoots(float a, float b, float c)
    {
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
        {
            t1 = -1;
            t2 = -1;
            return;
        }
        float sqrt_discr = sqrt( discr );
        t1 = (-b + sqrt_discr)/2.0/a;
        t2 = (-b - sqrt_discr)/2.0/a;
    }
    // Returns the lesser of the positive solutions, or a negative value if there was no positive solution.
    float getLesserPositive()
    {
        return ((0 < t1 && t1 < t2) || t2 < 0)?t1:t2;
    }
};


class Quadric: public Intersectable
{
public:
    float4x4 coeffs;
    
    Quadric(Material *material):
    Intersectable(material)
    {
        
        //initialize coeff matrix to make an ellipsoid
        coeffs = float4x4::identity();
        //make sure it will be visible (origin centered, modest radii)

    }
    void createSphere() {
        coeffs._33 = -1;
    }
    
    void createSmallYSphere() {
        coeffs._33 = -1;
        coeffs._11 = .25;
    }
    
    void createEllipsoid() {
        coeffs._33 = -1;
        coeffs._00 = 4;
    }
    
    Quadric* createCylinder() {

        
        coeffs._11 = 0;
        coeffs._33 = -1;
        return this;
    
    }
    
    Quadric* createWeirdCylinder() {
        
        
        coeffs._11 = 1.0;
       // coeffs._01 = .5;
      //  coeffs._22 = .5;
        coeffs._03 = .5;
        coeffs._33 = -1;
        return this;
        
    }
    Quadric* createCone() {
        coeffs._11 = -1;
        coeffs._33 =0;
        return this;
    }
    
    Quadric* createHyperboloid() {
        coeffs._11 = 1;
        coeffs._11 = -1;
        coeffs._22 = 1;
        coeffs._33= -1;

        return this;
    }
    
    
    Quadric* createHalfDome() {
        coeffs._11 = .5;
        coeffs._13 = 1;
        coeffs._33 = -1;

        return this;
    }
    
    Quadric* createRoundedBottom() {
        coeffs._11 = .5;
        coeffs._02 = 1;
        //OR for slightly less rounded do ._03
        coeffs._33 = -1;
        return this;
    }
    
    Quadric* createHorizontalCylinder() {
        coeffs._11 = .5;
        coeffs._00 = 0;
        coeffs._33 = -1;
        return this;
    }
    
    Quadric* transform(float4x4 t) {
        float4x4 tInverse = t.invert();
        coeffs = tInverse * coeffs * tInverse.transpose();
        return this;
    }
    
    
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        float4 rayDir4 = float4(ray.dir.x, ray.dir.y, ray.dir.z, 0);
        float term1 = rayDir4.dot(coeffs * rayDir4);
        float4 rayOrigin4 = float4(ray.origin.x, ray.origin.y, ray.origin.z, 1);
        float term2 = (rayDir4.dot(coeffs * rayOrigin4) + rayOrigin4.dot(coeffs * rayDir4));
        
        float term3 = rayOrigin4.dot(coeffs * rayOrigin4);
        
        return QuadraticRoots(term1, term2, term3);
        
    }
    float3 getNormalAt(float3 r)
    {
        float4 r4 = float4(r.x, r.y, r.z, 1);
        float4 term1 = coeffs * r4;
        float4 term2 = r4 * coeffs;
        float4 final = term1 + term2;
        
        float3 finalfinal = float3(final.x, final.y, final.z);
        return finalfinal.normalize();
    }
    
    Hit intersect(const Ray& ray)
    {
        //solve this with QuadraticRoots
        float t = solveQuadratic(ray).getLesserPositive();
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        return hit;
    }
    
    bool contains(float3 r)
    {
        float4 rhomo(r);
        // evaluate implicit eq
        float4 r4 = float4(r.x, r.y, r.z, 1);
        float result = r4.dot(coeffs * r4);
        // return true if negative
        if (result < 0)
            return true;
        else return false;
        // return false if positive
    }
    
    //infinite slab. ideal for clipping.
    Quadric* parallelPlanes() {
        coeffs = float4x4::identity();
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -1;
        return this;
    }
    
    Quadric* parallelPlanes2() {
        coeffs = float4x4::identity();
        coeffs._00 = 1;
        coeffs._11 = 0;
        coeffs._22 = -1;
        coeffs._33 = 0;
        return this;
    }
    
    
    //implement getNormalAt
    //r is a float3, but this works with float4s, so you have to append 1 for homogeneous then do product, then chop off fourth component, then normalize result
    
    //implement intersect--should be similar to Sphere intersect
    
};


//clipping: evaluate equation. if it's equal to 0, it's on the surface. if not, it's not on the surface.
class ClippedQuadric: public Intersectable
{
    Quadric shape = NULL;
    Quadric clipper = NULL;
 
public:
    ClippedQuadric(Material *material):
    Intersectable(material)
    
    {
    }
    
    float3 getNormalAt(float3 r)
    {
        float4 r4 = float4(r.x, r.y, r.z, 1);
        float4 term1 = shape.coeffs * r4;
        float4 term2 = r4 * shape.coeffs;
        float4 final = term1 + term2;
        
        float3 finalfinal = float3(final.x, final.y, final.z);
        return finalfinal.normalize();
    }
    
    ClippedQuadric* cylinder() {
        shape = *new Quadric(material);
        shape.createCylinder();
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* halfDome() {
        shape = *new Quadric(material);
        shape.createHalfDome();
        shape.transform(float4x4::rotation(float3(0.0,0.0,1.0),270));
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* weirdCylinder() {
        shape = *new Quadric(material);
        shape.createHyperboloid();
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* pawnCone() {
        shape = *new Quadric(material);
        shape.createCone();
        shape.transform(float4x4::translation(float3(0.0,1.0,0.0)) * float4x4::scaling(float3(.40,1.0,.40)));
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* bottomOfQueen() {
        shape = *new Quadric(material);
        shape.createHyperboloid();
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* oval() {
        shape = *new Quadric(material);
        shape.createSphere();
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    ClippedQuadric* topOfQueen() {
        shape = *new Quadric(material);
        shape.createEllipsoid();
        clipper = *new Quadric(material);
        clipper.parallelPlanes();
        return this;
    }
    
    
    ClippedQuadric* transform(float4x4 transformVector) {
        shape.transform(transformVector);
        clipper.transform(transformVector);
        return this;
    }

    
    
    
    //returns a bunch of data about the surface that the ray is intersecting. holds the position of the intersection, t (where along the ray it intersected), the normal, and the material.
    Hit intersect(const Ray& ray)
    {
        
        //solve quadratic equation for shape quadric
        QuadraticRoots twoRoots = shape.solveQuadratic(ray);
        float root1 = twoRoots.t1;
        float root2 = twoRoots.t2;
        
        //for both roots, check if intersection is within clipper quadric
        //use ray equation e + d dot t to see if intersection is within clipper quadric
        //substitute the coordinates of a point into the equation of the clipper
        
        float3 sol1 = ray.origin + ray.dir*root1;
        float3 sol2 = ray.origin + ray.dir*root2;
        
        //set root to negative value if not
        
        //do ray intersection for A
        //discard hits not in b
        if (!clipper.contains(sol1)) {
            twoRoots.t1 = -1;
        }
        if (!clipper.contains(sol2)) {
            twoRoots.t2 = -1;
        }
        
        float t = twoRoots.getLesserPositive();
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        return hit;
        
    }
};

class QueenCrown: public Intersectable
{
    Quadric shape = NULL;
    Quadric clipper1 = NULL;
    Quadric clipper2 = NULL;
    Quadric clipper3 = NULL;
    Quadric clipper4 = NULL;
public:
    QueenCrown(Material *material):
    Intersectable(material)
    
    {
    }
    
    float3 getNormalAt(float3 r)
    {
        float4 r4 = float4(r.x, r.y, r.z, 1);
        float4 term1 = shape.coeffs * r4;
        float4 term2 = r4 * shape.coeffs;
        float4 final = term1 + term2;
        
        float3 finalfinal = float3(final.x, final.y, final.z);
        return finalfinal.normalize();
    }
    
    QueenCrown* createCrown() {
        shape = *new Quadric(material);
        shape.createEllipsoid();
        
        clipper1 = *new Quadric(material);
        clipper1.parallelPlanes();
        clipper1.transform(float4x4::rotation(float3(0,0,1.0), -45));
        clipper1.transform(float4x4::scaling(float3(.25,.25,0.25)));
        clipper1.transform(float4x4::translation(float3(0.4, .4, 0)) );
        clipper2 = *new Quadric(material);
        clipper2.parallelPlanes();
        clipper2.transform(float4x4::rotation(float3(0,0,1.0), 45));
        clipper2.transform(float4x4::scaling(float3(.25,.25,.25)) );
        clipper2.transform(float4x4::translation(float3(-2.0, -2.0, 0)) );
        clipper3 = *new Quadric(material);
        clipper3.parallelPlanes();
        clipper3.transform(float4x4::translation(float3(0, 1.2, 0)) );
        clipper4 = *new Quadric(material);
        clipper4.createSphere();
        clipper4.transform(float4x4::scaling(float3(.35,.35,.35)) );
        clipper4.transform(float4x4::translation(float3(1.0, .9, 0)) );
        return this;
    }
    
    QueenCrown* transform(float4x4 transformVector) {
        shape.transform(transformVector);
        clipper1.transform(transformVector);
    
        clipper2.transform(transformVector);
        clipper3.transform(transformVector);
        clipper4.transform(transformVector);
        return this;
    }
    
    
    
    //returns a bunch of data about the surface that the ray is intersecting. holds the position of the intersection, t (where along the ray it intersected), the normal, and the material.
    Hit intersect(const Ray& ray)
    {
        
        //solve quadratic equation for shape quadric
        QuadraticRoots twoRoots = shape.solveQuadratic(ray);
        float root1 = twoRoots.t1;
        float root2 = twoRoots.t2;
        
        //for both roots, check if intersection is within clipper quadric
        //use ray equation e + d dot t to see if intersection is within clipper quadric
        //substitute the coordinates of a point into the equation of the clipper
        
        float3 sol1 = ray.origin + ray.dir*root1;
        float3 sol2 = ray.origin + ray.dir*root2;
        
        //set root to negative value if not
        
        //do ray intersection for A
        //discard hits not in b
        if (clipper1.contains(sol1)) {
            twoRoots.t1 = -1;
        }
        if (clipper1.contains(sol2)) {
            twoRoots.t2 = -1;
        }
        if (clipper2.contains(sol1)) {
            twoRoots.t1 = -1;
        }
        if (clipper2.contains(sol2)) {
            twoRoots.t2 = -1;
        }
        if (clipper3.contains(sol1)) {
            twoRoots.t1 = -1;
        }
        if (clipper3.contains(sol2)) {
            twoRoots.t2 = -1;
        }
        if (clipper4.contains(sol1)) {
            twoRoots.t1 = -1;
        }
        if (clipper4.contains(sol2)) {
            twoRoots.t2 = -1;
        }


        
        float t = twoRoots.getLesserPositive();
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        return hit;
        
    }
};


class Plane: public Intersectable
{
    float3 normal;
    float3 r0;
public:
    Plane(float3 normal, float3 r0, Material* material):
    Intersectable(material),
    normal(normal),
    r0(r0)
    {
        
    }
    
    Hit intersect(const Ray& ray)
    {
        
        float3 numerator = r0 - ray.origin;
        float finalNumerator = numerator.dot(normal);
        float denominator = ray.dir.dot(normal);
        
        float t= finalNumerator/denominator;
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = normal;
        return hit;
    }
};

//Plane class for chessboard
class RectangularPlane: public Intersectable
{
    float3 normal;
    float3 r0;
public:
    RectangularPlane(float3 normal, float3 r0, Material* material):
    Intersectable(material),
    normal(normal),
    r0(r0)
    {
        
    }
    
    Hit intersect(const Ray& ray)
    {
        
        float3 numerator = r0 - ray.origin;
        float finalNumerator = numerator.dot(normal);
        float denominator = ray.dir.dot(normal);
        
        float t= finalNumerator/denominator;

        
        //if position is outside bounds
        if (ray.origin.x + ray.dir.x * t > 4) {
            t = -1;
        }
        if (ray.origin.y + ray.dir.y * t >4) {
            t = -1;
        }
        if (ray.origin.z + ray.dir.z * t > 4) {
            t = -1;
        }
        if (ray.origin.x + ray.dir.x * t < -4) {
            t = -1;
        }
        if (ray.origin.y + ray.dir.y * t <-4) {
            t = -1;
        }
        if (ray.origin.z + ray.dir.z * t < -4) {
            t = -1;
        }
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = normal;
        return hit;
    }
};


class Scene
{
    Camera camera;

    DiffusePlane diffusePlaneMaterial;
    std::vector<Intersectable*> objects;
    HeadlightMaterial material;	// THIS NEEDS TO GO WHEN YOU USE A VECTOR OF MATERIALS
    
    std::vector<Material*> materials;
    std::vector<LightSource*> lightSources;
    HeadlightMaterial material2;
    Metal goldMetal;
    Metal silverMetal;
    PhongBlinn phongBlinnMaterialRed;
    PhongBlinn phongBlinnMaterialBlue;
    ProceduralTexturePhongBlinn proceduralTexturePhongBlinn;
    ProceduralTextureMetal proceduralTextureMetal;
    

    
public:
    Scene():
    diffusePlaneMaterial(),
    material2(float3(.8,.8,.8), float3(.8,.8,.8)),
    goldMetal(float3(.21,.485,1.29), float3(3.13, 2.23, 1.76)),
    silverMetal(float3(.15, .14, .13), float3(3.7, 3.11, 2.47)),
    phongBlinnMaterialRed(float3(1.0,1.0,1.0), float3(1.0,0.0,0.0), 15.0),
    phongBlinnMaterialBlue(float3(1.0,1.0,1.0), float3(0.0,0.0,1.0), 15.0),

    proceduralTexturePhongBlinn(float3(1.0,1.0,1.0), float3(0.0,0.0,1.0), 15.0),
    proceduralTextureMetal(float3(.21,.485,1.29), float3(3.13, 2.23, 1.76))

    {

        DirectionalLight *light1 = new DirectionalLight(float3(0.0, 1.0, 1.0), float3(1.0,1.0,1.0));
        lightSources.push_back(light1);
        
        //create chessboard
        RectangularPlane *rectPlane = new RectangularPlane(float3(0,.5,0), float3(0,-.9,0), &diffusePlaneMaterial);
       objects.push_back(rectPlane);

        //rook piece
//        ClippedQuadric *rook = new ClippedQuadric(&material);
//        rook->rook();
//        rook->transform(float4x4::scaling(float3(.20,.25,.25)) * float4x4::translation(float3(-0.75,-0.15,-0.2)));
        
        //pawn piece: sphere, oval, cone
        ClippedQuadric *firstPawn = new ClippedQuadric(&phongBlinnMaterialBlue);
        firstPawn->pawnCone();
        firstPawn->transform(float4x4::scaling(float3(.23,.28,.28)) * float4x4::translation(float3(1.0,-0.3,0.8)));
                             //* float4x4:: rotation(float3(1.0, 0.0, 0.0), M_PI/4));
        ClippedQuadric *pawnOval = new ClippedQuadric(&phongBlinnMaterialBlue);
        pawnOval->oval();
        pawnOval->transform(float4x4::scaling(float3(.15,.07,.25)) * float4x4::translation(float3(1.0,-0.10,0.8)));
        ClippedQuadric *pawnSphere = new ClippedQuadric(&phongBlinnMaterialBlue);
        pawnSphere->oval();
        pawnSphere->transform(float4x4::scaling(float3(.15,.15,.15)) * float4x4::translation(float3(1.0,0.10,0.8)));
        
        
        //knight piece
        ClippedQuadric *knightBottom = new ClippedQuadric(&proceduralTexturePhongBlinn);
        knightBottom->oval();
        knightBottom->transform(float4x4::scaling(float3(.24,.09,.21)) * float4x4::translation(float3(-1.40,-.85, .80)));
        ClippedQuadric *bottomKnightOval = new ClippedQuadric(&proceduralTexturePhongBlinn);
        bottomKnightOval->oval();
        bottomKnightOval->transform(float4x4::scaling(float3(.22,.10,.18)) * float4x4::translation(float3(-1.4,-.75,0.80)));
        ClippedQuadric *knightBodyBelly = new ClippedQuadric(&proceduralTexturePhongBlinn);
        ClippedQuadric *knightBody = new ClippedQuadric(&proceduralTexturePhongBlinn);
        knightBody->cylinder();
        knightBodyBelly->halfDome();
        knightBodyBelly->transform(float4x4::scaling(float3(.10,.12,.15)) * float4x4::translation(float3(-1.4,-.55,0.80)));
        knightBody->transform(float4x4::scaling(float3(.05,.12,.15)) * float4x4::translation(float3(-1.4,-.55,0.80)));
        //add half dome rotated 90
        ClippedQuadric *knightHead = new ClippedQuadric(&proceduralTexturePhongBlinn);
        knightHead->oval();
        knightHead->transform(float4x4::scaling(float3(.22,.16,.10)) * float4x4::translation(float3(-1.25,.90,.80)));
        knightHead->transform(float4x4::rotation(float3(0.0,0.0,1.0), 70));

        ClippedQuadric *test = new ClippedQuadric(&proceduralTexturePhongBlinn);
        test->weirdCylinder();
        test->transform(float4x4::rotation(float3(0,0,1.0), M_PI/2));

  
        //bishop piece
        ClippedQuadric *bishopBottom = new ClippedQuadric(&proceduralTextureMetal);
        bishopBottom->oval();
        bishopBottom->transform(float4x4::scaling(float3(.35,.30,.30)) * float4x4::translation(float3(-.5,-.95,-.5)));
        ClippedQuadric *bishopBody = new ClippedQuadric(&proceduralTextureMetal);
        bishopBody->bottomOfQueen();
        bishopBody->transform(float4x4::scaling(float3(.15,.40,.15)) * float4x4::translation(float3(-.5,-.45,-.5)));
        ClippedQuadric *bishopBigOval = new ClippedQuadric(&proceduralTextureMetal);
        bishopBigOval->oval();
        bishopBigOval->transform(float4x4::scaling(float3(.30,.12,.25)) * float4x4::translation(float3(-0.5,-.11,-.5)));
        ClippedQuadric *bishopSmallOval1 = new ClippedQuadric(&proceduralTextureMetal);
        bishopSmallOval1->oval();
        bishopSmallOval1->transform(float4x4::scaling(float3(.17,.05,.15)) * float4x4::translation(float3(-0.5,.02,-.5)));
        ClippedQuadric *bishopSmallOval2 = new ClippedQuadric(&proceduralTextureMetal);
        bishopSmallOval2->oval();
        bishopSmallOval2->transform(float4x4::scaling(float3(.17,.05,.15)) * float4x4::translation(float3(-0.5,0.1,-.5)));
        ClippedQuadric *bishopHead = new ClippedQuadric(&proceduralTextureMetal);
        bishopHead->oval();
        bishopHead->transform(float4x4::scaling(float3(.20,.30,.15)) * float4x4::translation(float3(-0.5,0.25,-.5)));
        ClippedQuadric *bishopTinyTop = new ClippedQuadric(&proceduralTextureMetal);
        bishopTinyTop->oval();
        bishopTinyTop->transform(float4x4::scaling(float3(.05,.05,.05)) * float4x4::translation(float3(-0.5,0.60,-.5)));
        
        //queen piece
        ClippedQuadric *bottomOfQueen = new ClippedQuadric(&goldMetal);
        bottomOfQueen->bottomOfQueen();
        bottomOfQueen->transform(float4x4::scaling(float3(.18,.40,.25)) * float4x4::translation(float3(0.35,-.35,0.9)));
//        ClippedQuadric *topOfQueen = new ClippedQuadric(&goldMetal);
//        topOfQueen->topOfQueen();
//        topOfQueen->transform(float4x4::scaling(float3(.40,.30,.25)) * float4x4::translation(float3(0.35,0.05,0.9)));
        QueenCrown *queenCrown = new QueenCrown(&goldMetal);
        queenCrown->createCrown();
        queenCrown->transform(float4x4::scaling(float3(.50,.50,.25))* float4x4::translation(float3(0.45, .55, 0)));
        
        
//        ClippedQuadric *queenOval = new ClippedQuadric(&metalMaterial);
//        queenOval->oval();
//        queenOval->transform(float4x4::scaling(float3(.24,.10,.25)) * float4x4::translation(float3(0.35,-.2,0.9)));
        
 //       objects.push_back(rook);
        objects.push_back(firstPawn);
        objects.push_back(bottomOfQueen);
   //     objects.push_back(topOfQueen);
        objects.push_back(queenCrown);
        objects.push_back(pawnOval);
        objects.push_back(pawnSphere);
        objects.push_back(knightBottom);
        objects.push_back(knightHead);
        objects.push_back(knightBody);
        objects.push_back(knightBodyBelly);
        objects.push_back(bottomKnightOval);
        objects.push_back(bishopBottom);
        objects.push_back(bishopBody);
        objects.push_back(bishopBigOval);
        objects.push_back(bishopSmallOval1);
        objects.push_back(bishopSmallOval2);
        objects.push_back(bishopHead);
        objects.push_back(bishopTinyTop);
       //objects.push_back(test);

    }
    ~Scene()
    {
        // UNCOMMENT THESE WHEN APPROPRIATE
        //for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
        //	delete *iMaterial;
        //for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
        //	delete *iObject;
    }
    
public:
    Camera& getCamera()
    {
        return camera;
    }
    
    Hit firstIntersect(Ray ray) {
        Hit bestHit;
        bestHit.t = FLT_MAX;
        for(Intersectable* obj : objects)
        {
            Hit hit = obj->intersect(ray);
            if(hit.t > 0 && hit.t < bestHit.t) {
                bestHit = hit;
            }
        }
        if (bestHit.t == FLT_MAX) {
            bestHit = Hit();
        }
        return bestHit;
    }
    
    float3 trace(const Ray& ray, int depth)

    {
        if (depth > 1000) {
            return float3(0,0,0);
        }
        Hit hit = firstIntersect(ray);
        if(hit.t < 0)
            return float3(1, 1, 1);
        
        float3 sum = float3(0,0,0);
        Metal* metal = dynamic_cast< Metal*>(hit.material);
        ProceduralTextureMetal *proceduralTextureMetal =dynamic_cast< ProceduralTextureMetal*>(hit.material);
        //TODO: add switch statement of metal
        if (metal != NULL) {
            Metal::Event e = metal->evaluateEvent(ray.dir, hit.normal);
            float3 reflectedRay = e.reflectionDir;
            float3 reflectance = e.reflectance;
            //if hit normal is pointing to other side
            //dot product of normal and viewDir
            float3 normal = hit.normal;
            if (ray.dir.dot(normal) < 0) {
                
            }
            else {
                normal = -normal;
            }
            
            Ray reflectionRay(hit.position + normal *.01, reflectedRay);
            return trace(reflectionRay, depth+1) * reflectance;
            
        }
        else if (proceduralTextureMetal != NULL) {
            ProceduralTextureMetal::Event e = proceduralTextureMetal->evaluateEvent(ray.dir, hit.normal, hit.position);
            float3 reflectedRay = e.reflectionDir;
            float3 reflectance = e.reflectance;
            //if hit normal is pointing to other side
            //dot product of normal and viewDir
            float3 normal = hit.normal;
            if (ray.dir.dot(normal) < 0) {
                
            }
            else {
                normal = -normal;
            }
            
            Ray reflectionRay(hit.position + normal *.01, reflectedRay);
            return trace(reflectionRay, depth+1) * reflectance;

        }
        
        else {
            for (int i =0; i < lightSources.size(); i++) {
                
                float3 currentLightDir = lightSources.at(i)->getLightDirAt(hit.position);
                
                Ray shadowRay(hit.position + hit.normal* .01,currentLightDir);
                
                float3 currentPowerDensity = lightSources.at(i)->getPowerDensityAt(hit.position);
                Hit shadowHit = firstIntersect(shadowRay);
                if (shadowHit.t < 0 || shadowHit.t > lightSources.at(i)->getDistanceFrom(hit.position)) {
                    sum += hit.material->shade(hit.normal, -ray.dir, currentLightDir, currentPowerDensity, hit.position);
                    
                }
                
            }
            //if it exceeds depth value, just return black
        }
        
        
        return sum;
        //return hit.material->getColor(hit.position, hit.normal, -ray.dir);
        
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
    static unsigned int iPart = 0;
    
    if(iPart >= 64)
        return false;
    for(int j = iPart; j < screenHeight; j+=64)
    {
        for(int i = 0; i < screenWidth; i++)
        {
            float3 pixelColor = float3(0, 0, 0);
            float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );
            
            Camera& camera = scene.getCamera();
            Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
            
            image[j*screenWidth + i] = scene.trace(ray, 0);
        }
    }
    iPart++;
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen
    
    if(computeImage())
        glutPostRedisplay();
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
    
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(screenWidth, screenHeight);				// startup window size
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
    
    glutCreateWindow("Ray caster");				// application window is created and displayed
    
    glViewport(0, 0, screenWidth, screenHeight);
    
    glutDisplayFunc(onDisplay);					// register callback
    
    glutMainLoop();								// launch event handling loop
    
    return 0;
}

