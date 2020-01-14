Shader "Custom/RT"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma enable_d3d11_debug_symbols

            #include "UnityCG.cginc"

            static const int lightCount = 2;
            static const int sphereCount = 4;
            static const int planeCount = 1;
            static const int cylinderCount = 2;
            static const int maxBounce = 10;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };


            struct PointLight{
                float3 position;
                float3 color;
            };

            struct Material{
                float3 diffuse;
                float3 specular;
                float glossiness;
                // Reflectance
                float kr;
                // Transmittance
                float kt;
                float IOR;
            };

            struct Sphere{
                float3 position;
                float radius;
                Material material;
                bool enable;
            };

            struct Plane{
                float3 normal;
                float d;
                Material material;
                bool enable;
            };

            struct Cylinder{
                float3 position;
                float3 direction;
                float radius;
                Material material;
                bool enable;
            };

            struct Scene{
                float3 ambient;
                PointLight lights[lightCount];
                Sphere spheres[sphereCount];
                Plane planes[planeCount];
                Cylinder cylinders[cylinderCount];
            };

            struct Ray{
                float3 origin;
                float3 direction;
            };

            // Contains all information pertaining to a ray/object intersection
            struct HitInfo{
                bool hit;
                bool inside;
                float t;
                float3 position;
                float3 normal;
                Material material;
            };

            HitInfo getEmptyHit(){
                Material emptyMat;
                emptyMat.diffuse = float3(0, 0, 0);
                emptyMat.specular = float3(0, 0, 0);
                emptyMat.glossiness = 0;
                emptyMat.kr = 0;
                emptyMat.kt = 0;
                emptyMat.IOR = 1;

                HitInfo empHit;
                empHit.hit = false;
                empHit.inside = false;
                empHit.t = 0;
                empHit.position = float3(0, 0, 0);
                empHit.normal = float3(0, 0, 0);
                empHit.material = emptyMat;

                return empHit;
            }

            // Sorts the two t values such that t1 is smaller than t2
            void sortT(inout float t1,inout float t2){
                // Make t1 the smaller t
                if(t2<t1){
                    float temp=t1;
                    t1=t2;
                    t2=temp;
                }
            }

            // Tests if t is in an interval
            bool isTInInterval(float t, float tMin, float tMax){
                return t>tMin&&t<tMax;
            }

            // Get the smallest t in an interval
            bool getSmallestTInInterval(float t0,float t1,float tMin, float tMax, inout float smallestTInInterval){
                sortT(t0,t1);
                // As t0 is smaller, test this first
                if(isTInInterval(t0,tMin,tMax)){
                    smallestTInInterval=t0;
                    return true;
                }
                
                // If t0 was not in the interval, still t1 could be
                if(isTInInterval(t1,tMin,tMax)){
                    smallestTInInterval=t1;
                    return true;
                }
                // none was
                return false;
            }

            // Get the smallest t in an interval
            bool getBiggestTInInterval(float t0,float t1, float tMin, float tMax,inout float biggestTInInterval){
                sortT(t0,t1);
                // As t0 is smaller, test this first
                if(isTInInterval(t1,tMin,tMax)){
                    biggestTInInterval=t1;
                    return true;
                }
                
                // If t0 was not in the interval, still t1 could be
                if(isTInInterval(t0,tMin,tMax)){
                    biggestTInInterval=t0;
                    return true;
                }
                // none was
                return false;
            }

            HitInfo intersectSphere(Ray ray, Sphere sphere, float tMin, float tMax){
                
                float3 to_sphere=ray.origin-sphere.position;
                
                float a=dot(ray.direction,ray.direction);
                float b=2.*dot(ray.direction,to_sphere);
                float c=dot(to_sphere,to_sphere)-sphere.radius*sphere.radius;
                float D=b*b-4.*a*c;
                if(D>0.)
                {
                    float t0=(-b-sqrt(D))/(2.*a);
                    float t1=(-b+sqrt(D))/(2.*a);
                    
                    float smallestTInInterval;
                    if(!getSmallestTInInterval(t0,t1,tMin,tMax,smallestTInInterval)){
                        return getEmptyHit();
                    }
                    HitInfo hitInfo;

                    hitInfo.position =ray.origin+smallestTInInterval*ray.direction;
                    if (length(ray.origin - sphere.position) < sphere.radius + .0001)
                    {
                        hitInfo.normal = -normalize(hitInfo.position - sphere.position);
                        hitInfo.inside = true;
                    }
                    else {
                        hitInfo.normal = normalize(hitInfo.position - sphere.position);
                        hitInfo.inside = false;
                    }
                    
                    hitInfo.hit = true;
                    hitInfo.t = smallestTInInterval;
                    hitInfo.material = sphere.material;
                    return hitInfo;
                }
                return getEmptyHit();
            }

            HitInfo intersectPlane(const Ray ray,const Plane plane,const float tMin,const float tMax){
                // Add your plane intersection code here
                //Let the intersection point be P(t) and a point on the plane be V0
                //P(t) - V0 = ray.direction * t - (V0 - ray.origin)
                //Since P(t) - V0 is a line on the Plance, dot(P(t) - V0, plane.normal) == 0
                // dot(ray.direction * t - (V0 - ray.origin), plane.normal) == 0
                // t = dot(V0 - ray.origin, plane.normal) / dot(ray.direction, plane.normal)
                float norDotdir=dot(plane.normal,ray.direction);
                float3 normal;
                float t;
                
                //if dot product is -ve, the ray hits from the front side of the plane, which is good
                if(norDotdir<0.){
                    normal=plane.normal;
                }
                //if dot product is +ve, the ray hits from the back side of the plane, so we need to flip the plane by reversing the normal
                else if(norDotdir>0.){
                    normal=-plane.normal;
                    norDotdir=-norDotdir;
                }
                //if dot product is 0, the ray is parallel to the plane, which would be considered not hitted.
                else
                    return getEmptyHit();
                //Use the above equation to calculate t
                t=(dot(normal,float3(0.,-plane.d,0.)-ray.origin))/norDotdir;
                //Make sure the t is within range
                if(!isTInInterval(t,tMin,tMax))
                return getEmptyHit();
                //Calculate the hit position with the line equation
                float3 hitPosition=ray.origin+t*ray.direction;
                HitInfo hitInfo;
                hitInfo.hit = true;
                hitInfo.inside = false;
                hitInfo.t = t;
                hitInfo.position = hitPosition;
                hitInfo.normal = normal;
                hitInfo.material = plane.material;
                return hitInfo;
            }

            float lengthSquared(float3 x){
                return dot(x,x);
            }

            HitInfo intersectCylinder(const Ray ray,const Cylinder cylinder,const float tMin,const float tMax){
                // Add your cylinder intersection code here
                // Let the equation for a line be R(t) p0 + v0t, where p0 is a point on the line and v0 be the directional vector
                // and the equation for a cylinder be (q - p1 - dot(v1, (q - p1)) * v1)^2 = r^2, where q is a point on the cylinder surface,
                // v1 and p1 is the direction and point representing the line which the cylinder is oriented,
                // and r is the radius of the cylinder
                // Solving the equations gives: (v0-dot(v1, v0) * v1)^2 * t^2 + 2*dot(v0-dot(v1, v0) * v1, p0-p1 - dot(p0-p1, v1) * v1) * t + (p0-p1 - dot(p0-p1, v1) * v1)^2 - r^2,
                // which is just a quadratic equation of t
                float3 RC=(ray.origin-cylinder.position);
                float3 normalVec=cross(ray.direction,cylinder.direction);
                float ln=length(normalVec);
                
                //if the ray and the cylinder's direction is nearly parallel, ignore the hit
                if(ln<.01)
                {
                    return getEmptyHit();
                }
                float a=lengthSquared(ray.direction-dot(ray.direction,cylinder.direction)*cylinder.direction);
                float b=2.*dot(ray.direction-dot(ray.direction,cylinder.direction)*cylinder.direction,RC-dot(RC,cylinder.direction)*cylinder.direction);
                float c=lengthSquared(RC-dot(RC,cylinder.direction)*cylinder.direction)-cylinder.radius*cylinder.radius;
                float D=b*b-4.*a*c;
                if (D > 0.)
                {
                    float t0 = (-b - sqrt(D)) / (2. * a);
                    float t1 = (-b + sqrt(D)) / (2. * a);

                    float smallestTInInterval;
                    //Get the closer t to the origin in the +ve direction
                    if (!getSmallestTInInterval(t0, t1, tMin, tMax, smallestTInInterval)) {
                        return getEmptyHit();
                    }
                    float3 hitPosition = ray.origin + smallestTInInterval * ray.direction;
                    // In order to figure out the normal, we need to find a point on the cylinder central line, p2,
                    // such that (hitPosition - p2) is perpendicular to the cylinder direction
                    // We do so by projecting the vector (hitPosition - cylinder.position) onto the cylinder direction vector
                    float3 p2 = cylinder.position + cylinder.direction * dot((hitPosition - cylinder.position), cylinder.direction);
                    // The normal is simply the normalized (hitPosition - p2)
                    float3 normal = normalize(hitPosition - p2);
                    HitInfo hitInfo;

                    // If the ray direction and the normal is pointing to the kind of same direction,
                    // i.e. the angle between them is smaller than 90 degree,
                    // the ray hits from inside the cylinder
                    if (dot(normal, ray.direction) > 0.) {
                        //flip the normal
                        normal = -normal;
                        hitInfo.inside = true;
                    }
                    else
                        hitInfo.inside = false;

                    hitInfo.hit = true;
                    hitInfo.t = smallestTInInterval;
                    hitInfo.position = hitPosition;
                    hitInfo.normal = normal;
                    hitInfo.material = cylinder.material;
                    return hitInfo;
                }
                else
                    return getEmptyHit();
            }

            HitInfo getBetterHitInfo(const HitInfo oldHitInfo,const HitInfo newHitInfo){
                if(newHitInfo.hit)
                if(newHitInfo.t<oldHitInfo.t)// no need to test for the interval, this has to be done per-primitive
                return newHitInfo;
                return oldHitInfo;
            }

            HitInfo intersectScene(const Scene scene,const Ray ray,const float tMin,const float tMax){
                HitInfo bestHitInfo;
                bestHitInfo.t=tMax;
                bestHitInfo.hit=false;
                for (int i = 0; i < cylinderCount; ++i) {
                    if(scene.cylinders[i].enable)
                        bestHitInfo = getBetterHitInfo(bestHitInfo, intersectCylinder(ray, scene.cylinders[i], tMin, tMax));
                }
                for (int i = 0; i < sphereCount; ++i) {
                    if (scene.spheres[i].enable)
                        bestHitInfo = getBetterHitInfo(bestHitInfo, intersectSphere(ray, scene.spheres[i], tMin, tMax));
                }
                for (int i = 0; i < planeCount; ++i) {
                    if (scene.planes[i].enable)
                        bestHitInfo = getBetterHitInfo(bestHitInfo, intersectPlane(ray, scene.planes[i], tMin, tMax));
                }
                
                return bestHitInfo;
            }

            float3 shadeFromLight(
                const Scene scene,
                const Ray ray,
                const HitInfo hit_info,
            const PointLight light)
            {
                float3 hitToLight=light.position-hit_info.position;
                
                float3 lightDirection=normalize(hitToLight);
                float3 viewDirection=normalize(hit_info.position-ray.origin);
                float3 reflectedDirection=reflect(viewDirection,hit_info.normal);
                float diffuse_term=max(0.,dot(lightDirection,hit_info.normal));
                float specular_term=pow(max(0.,dot(lightDirection,reflectedDirection)),hit_info.material.glossiness);
                
                // Put your shadow test here
                float visibility=1.;
                // Create a ray from the object hit point to the light.
                Ray shadowRay;
                shadowRay.origin = hit_info.position;
                shadowRay.direction = lightDirection;
                // Notice that I use the distance beteen the object hit point and the light for the maxium length of the intersection check
                // If I use a length longer that, it could return a hit with an object behind the light
                // Also, if I choose a minimum length of 0, shaow-acne would occur, where dots appear on the surfaces the of objects.
                // It is caused by numerical precision errors which make the ray intersect the surface from which it is cast.
                HitInfo shadowHitInfo=intersectScene(scene,shadowRay,.001,length(hitToLight));
                // If it return a hit, that means something is blocking the light
                /* if(shadowHitInfo.hit){
                    // A simple shadow can end here and just return zero, so that no color is returned from this light
                    visibility = 0.0;
                    // But some object is transparent, and light can pass through them
                    if(shadowHitInfo.material.kt>0.0){
                        // The more perpendicular between the normal and the light direction, the higher the visibility
                        float x = abs(dot(shadowRay.direction, shadowHitInfo.normal)) * shadowHitInfo.material.kt;
                        // Exaggerated the maximum visibility to fake the converging light, only works for convex shape
                        visibility = clamp(pow(x, 6.0) * 7.0, 0.0, 2.5);
                    }
                    // This isn't near perfect though. What if there is a third object blocking the light from both the original object and
                    // the second transparent object. I think I need to keep casting shadow rays until nothing is hitted between the light
                    // and the n-th object, and then calculate the remained light received after passing through n objects.
                } */
                for (int i=0;i < 3;i++){
                    if(shadowHitInfo.hit){
                        // But some object is transparent, and light can pass through them
                        if(shadowHitInfo.material.kt>0.0){
                            // The more perpendicular between the normal and the light direction, the higher the visibility
                            float x = abs(dot(shadowRay.direction, shadowHitInfo.normal)) * shadowHitInfo.material.kt;
                            // Exaggerated the maximum visibility to fake the converging light, only works for convex shape
                            //visibility *= clamp(pow(x, 6.0) * 7.0, 0.0, 2.7);
                            visibility *= clamp(pow(x, 6.0) * 7.0, 0.0, 1.);
                            shadowRay.origin = shadowHitInfo.position;
                            shadowHitInfo = intersectScene(scene, shadowRay, 0.001, length(light.position - shadowHitInfo.position));
                            shadowRay.origin = shadowHitInfo.position;
                            shadowHitInfo = intersectScene(scene, shadowRay, 0.001, length(light.position - shadowHitInfo.position));
                        }
                        else{
                            visibility = 0.0;
                            break;
                        }
                    }
                }
                
                Ray mirrorRay;
                mirrorRay.origin = hit_info.position;
                mirrorRay.direction = reflect(lightDirection, hit_info.normal);
                HitInfo mirrorHitInfo = intersectScene(scene, mirrorRay, 0.001, 100000.0);
                
                return visibility *
                light.color * (
                    specular_term * hit_info.material.specular
                    +
                    diffuse_term * hit_info.material.diffuse
                );
            }

            float3 background(const Ray ray) {
                // A simple implicit sky that can be used for the background
                //return float3(ray.direction.y, 0, 0);
                //return float3(0.0) + float3(0.0, 1., 0.0) * max(0.0, ray.direction.y); ///2.0);
                return float3(0.2,.2,.2) + float3(0.8, 0.6, 0.5) * max(0.0, ray.direction.y);
            }

            // It seems to be a WebGL issue that the third parameter needs to be inout instead of const on Tobias' machine
            float3 shade(const Scene scene, const Ray ray, inout HitInfo hitInfo) {
                
                if(!hitInfo.hit) {
                    //return float3(1);
                    return background(ray);
                }
                float3 shading = scene.ambient * hitInfo.material.diffuse;
                for (int i = 0; i < lightCount; ++i) {
                    shading += shadeFromLight(scene, ray, hitInfo, scene.lights[i]);
                }
                return shading;
            }

            Ray getFragCoordRay(const float2 frag_coord) {
                float sensorDistance = 1.0;
                float2 sensorMin = float2(-1, 0.5);
                float2 sensorMax = float2(1, -0.5);
                float2 pixelSize = (sensorMax - sensorMin) / float2(800, 400);
                float3 origin = float3(0, 0, sensorDistance);
                float3 direction = normalize(float3(sensorMin + pixelSize * frag_coord, -sensorDistance));
                Ray ray;
                ray.origin = origin;
                ray.direction = direction;
                return ray;
            }

            // Added field for the ior of the object 
            float fresnel(const float3 viewDirection, const float3 normal, float n2) {
                // Put your code to compute the Fresnel effect here
                float n1 = 1.;
                float cosi = abs(dot(viewDirection, normal));
                
                // From Schlick's approximation on wikipedia
                // R = R0 + (1 - R0)*(1-cosi)^5, where R0 = ((n1 - n2)/(n1+n2))^2
                /* float R0 = ((n1 - n2)/(n1+n2)) * ((n1 - n2)/(n1+n2));
                return R0 + (1. - R0)*pow(1.-cosi, 5.); */
                // Actually I can't notice any difference between the Schlick's approximation and the Fresnel equations below...
                
                // From wikipedia, the effective reflectivity of unpolarized light can be calculated by:
                // R = 1/2 * (Rs + Rp), where Rs and Rp is the reflectance for s-polarized light and p-polarized light
                // and Rs = ((n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost))^2
                // and Rp = ((n1 * cost - n2 * cosi) / (n1 * cost + n2 * cosi))^2
                
                // As both vectors are normalized, I won't bother to divide it by their lengths
                float sini = sqrt(max(0.0, 1.0 - cosi * cosi));
                float sint = n1 / n2 * sini;
                
                // Normal incidence, where there is no distinction between s and p polarization
                if(sini == 0. && sint == 0.){
                    return ((n1 - n2) / (n1+ n2)) * ((n1 - n2) / (n1+ n2));
                }
                // If the refracted angle is larger than 90 degree, total internal reflection has occurred
                // If the angle is 90, also return 1.0 to avoid trouble
                else if (sint >= 1.0) {
                    return 1.0;
                }
                else {
                    // Brewster's angle
                    float sinb = sin(atan(n2/n1));
                    float cost = sqrt(max(0.0, 1.0 - sint * sint));
                    float Rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost);
                    Rs *= Rs;
                    float Rp;
                    // If angle of incidence is at the Brewster's angle, Rp goes to zero
                    if (abs(sinb - sini) == 0.0)
                    Rp = 0.;
                    else{
                        Rp = (n1 * cost - n2 * cosi) / (n1 * cost + n2 * cosi);
                        Rp *= Rp;
                    }
                    
                    return (Rs + Rp) / 2.0;
                }
            }

            float3 bounce(const Scene scene, const Ray currentRay, const HitInfo currentHitInfo, float weight, int step) {
                float3 result = float3(0, 0, 0);
                //Determine whether the ray is inside the hitted object
                float IORi;
                float IORt;
                //If the ray is inside the object, I assume the ray is refracting out to air
                if (currentHitInfo.inside) {
                    IORi = currentHitInfo.material.IOR;
                    IORt = 1.0;
                }
                else {
                    IORi = 1.0;
                    IORt = currentHitInfo.material.IOR;
                }
                float eta = IORi / IORt;
                // Fresnel
                float f = fresnel(currentRay.direction, currentHitInfo.normal, eta);

                if (f > 0.05) {
                    // Compute the reflection
                    float reflectionWeight = weight * currentHitInfo.material.kr * f;

                    Ray nextRay;
                    nextRay.origin = currentHitInfo.position;
                    nextRay.direction = currentRay.direction - 2.0 * dot(currentHitInfo.normal, currentRay.direction) * currentHitInfo.normal;
                    HitInfo nextHitInfo = intersectScene(scene, nextRay, 0.001, 10000.0);
                    if (nextHitInfo.hit)
                    {
                        result += reflectionWeight * shade(scene, nextRay, nextHitInfo);
                        if (step < maxBounce && reflectionWeight > 0.005)
                            result += bounce(scene, nextRay, nextHitInfo, reflectionWeight, ++step);
                    }
                }
                if (1 - f > 0.05) {
                    // Compute the refraction
                    float refractionWeight = weight * currentHitInfo.material.kt * (1-f);
                    
                    Ray nextRay;
                    float cosa = dot(-currentRay.direction, currentHitInfo.normal);
                    float root = 1.0 + eta * eta * (cosa * cosa - 1.0);
                    float w;
                    if (root < 0.0) {
                    //    // total internal reflection
                    //    nextRay.origin = currentHitInfo.position;
                    //    nextRay.direction = reflect(currentRay.direction, currentHitInfo.normal);
                    //    w = reflectionWeight;
                    }
                    else {
                        // refraction
                        nextRay.origin = currentHitInfo.position;
                        nextRay.direction = -eta * -currentRay.direction + currentHitInfo.normal * (eta * cosa - sqrt(root));
                        //refract(currentRay.direction, currentHitInfo.normal, eta));
                        //w = refractionWeight;
                        //refractionWeight = w;

                        HitInfo nextHitInfo = intersectScene(scene, nextRay, .001, 10000.);
                        if (nextHitInfo.hit)
                        {
                            result += refractionWeight * shade(scene, nextRay, nextHitInfo);
                            if (step < maxBounce && refractionWeight > 0.005)
                                result += bounce(scene, nextRay, nextHitInfo, refractionWeight, ++step);
                        }
                    }
                }
                return result;
            }

            float3 colorForFragment(const Scene scene, const float2 fragCoord) {
                
                Ray initialRay = getFragCoordRay(fragCoord);
                
                HitInfo initialHitInfo = intersectScene(scene, initialRay, 0.001, 10000.0);
                float3 result = shade(scene, initialRay, initialHitInfo);
                if (initialHitInfo.hit) {
                    //int step = 0;
                    //while (step < MaxBounce) {

                    //    step++;
                    //}
                    result += bounce(scene, initialRay, initialHitInfo, 1., 0);
                }
                
                return result;
            }
                
            // --------------------------------------------Table of Material----------------------------------------------------
            // |  Material      |       Diffuse          | Specular    |   Glossiness   | Reflectance | Transmittance |   IOR  |
            // -----------------------------------------------------------------------------------------------------------------
            // |     Paper      |       White(1)         |   None(0)   |    Low, 1      |  Low, 0     |   Low, 0.3    |  1.557 |
            // -----------------------------------------------------------------------------------------------------------------
            // |                | Paper is white (normally) and has scattered reflection, hence the low glosinness. A thin     |
            // |      Why       | sheet of paper can transmitt a small amount of light. The IOR is collected from a paper      |
            // |                | using streakmem measurement.                                                                 |
            // -----------------------------------------------------------------------------------------------------------------
            // |  Steel Mirror  |     Dark Gray(0.1)     |   High(1)   | very High, 100 |  High, 1    |   None, 0     |   2.5  |
            // -----------------------------------------------------------------------------------------------------------------
            // |      Why       | Steel Mirror is Gray (normally) and highly reflective with its very smooth surface. It has   |
            // |                | negligible transmittance and high IOR.                                                       |                                                      |
            // -----------------------------------------------------------------------------------------------------------------
            // |     Glass      |        None(0)         |   High(1)   | very High, 100 |  High, 1    |   High, 1     |   1.5  |
            // -----------------------------------------------------------------------------------------------------------------
            // |      Why       | Glass is fully transparent (I assumed) and highly reflective with its very smooth surface.   |
            // |                | It also has high transmittance and a comparatively low IOR.                                  |
            // -----------------------------------------------------------------------------------------------------------------
            // | Yellow plastic | Yellow-ish(1, 0.25, 0) |   High(1)   |     High, 10   |  High, 1    |   None, 0     |  1.46  |
            // -----------------------------------------------------------------------------------------------------------------
            // |                | Yellow plastic is opaque (I assumed) and highly reflective with its smooth surface (but not  |
            // |      WHy       | as smooth as the glass or steel). I think its color is more close to the orange color, so I  |
            // |                | set the diffuse with higher red value. It also has a low IOR.                                |
            // -----------------------------------------------------------------------------------------------------------------
                
            Material getDefaultMaterial(){
                // Update the default material call to match the new parameters of Material
                Material mat;
                mat.diffuse = float3(.3, .3, .3);
                mat.specular = float3(0, 0, 0);
                mat.glossiness = 1;
                mat.kr = 0;
                mat.kt = 0;
                mat.IOR = 0;
                return mat;
            }
                
            Material getPaperMaterial(){
                // Replace by your definition of a paper material
                Material mat;
                mat.diffuse = float3(1, 1, 1);
                mat.specular = float3(0, 0, 0);
                mat.glossiness = .3;
                mat.kr = 0;
                mat.kt = .3;
                mat.IOR = 1.557;
                return mat;
            }
                
            Material getPlasticMaterial(){
                // Replace by your definition of a plastic material
                Material mat;
                mat.diffuse = float3(.25, 1., .2);
                mat.specular = float3(1,1,1);
                mat.glossiness = 10;
                mat.kr = 1;
                mat.kt = 0;
                mat.IOR = 1.46;
                return mat;
            }
                
            Material getGlassMaterial(){
                // Replace by your definition of a glass material
                Material mat;
                mat.diffuse = float3(0, 0, 0);
                mat.specular = float3(1, 1, 1);
                mat.glossiness = 100;
                mat.kr = 1;
                mat.kt = 1;
                mat.IOR = 1.5;
                return mat;
            }
                
            Material getSteelMirrorMaterial(){
                // Replace by your definition of a steel mirror material
                Material mat;
                mat.diffuse = float3(.1,.1,.1);
                mat.specular = float3(1, 1, 1);
                mat.glossiness = 100;
                mat.kr = 1;
                mat.kt = 0;
                mat.IOR = 2.5;
                return mat;
            }

            Material getPlasticMirrorMaterial() {
                // Replace by your definition of a steel mirror material
                Material mat;
                mat.diffuse = float3(1, .1, 0);
                mat.specular = float3(1, 1, 1);
                mat.glossiness = 100;
                mat.kr = 1;
                mat.kt = 0;
                mat.IOR = 2.5;
                return mat;
            }
                
            float3 tonemap(const float3 radiance){
                const float invGamma = 0.5;
                return pow(radiance,float3(invGamma, invGamma, invGamma));
            }
                
            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                //o.vertex = v.vertex;
                o.uv = v.uv;
                return o;
            }

            sampler2D _MainTex;

            fixed4 frag (v2f i) : SV_Target
            {
                // Setup scene
                Scene scene;
                scene.ambient = float3(.2,.15,.2);

                // Lights
                scene.lights[0].position = float3(5,15,-5);
                scene.lights[0].color = .5 * float3(.9,.5,.1);

                scene.lights[1].position = float3(-15,5,2);
                scene.lights[1].color = .5 * float3(.1,.3,1.);

                // Primitives
                scene.spheres[0].position = float3(10,-5,-16);
                scene.spheres[0].radius = 6.;
                scene.spheres[0].material = getPaperMaterial();
                scene.spheres[0].enable = true;
                //getGlassMaterial();

                scene.spheres[1].position = float3(-7,-1,-13);
                scene.spheres[1].radius = 4.;
                scene.spheres[1].material = getPlasticMaterial();
                scene.spheres[1].enable = true;
                //getGlassMaterial();

                scene.spheres[2].position = float3(5 * _CosTime.w ,.5, -10 + 5 * _SinTime.w); //float3(0 + _SinTime.w,.5,-5)
                scene.spheres[2].radius = 2.;
                scene.spheres[2].material = getGlassMaterial();
                scene.spheres[2].enable = true;

                scene.spheres[3].position = float3(5 * cos(_Time.y+90), .5, -10 + 5 * sin(_Time.y + 90));
                scene.spheres[3].radius = 2.;
                scene.spheres[3].material = getPlasticMirrorMaterial();
                scene.spheres[3].enable = true;

                scene.planes[0].normal = float3(0,1,0);
                scene.planes[0].d = 4.5;
                scene.planes[0].material = getSteelMirrorMaterial();
                scene.planes[0].enable = true;
                //getGlassMaterial();

                scene.cylinders[0].position = float3(-1,1,-18);
                scene.cylinders[0].direction = normalize(float3(-1,2,-1));
                scene.cylinders[0].radius = 1.5;
                scene.cylinders[0].material = getPaperMaterial();
                scene.cylinders[0].enable = true;
                //getGlassMaterial();

                scene.cylinders[1].position = float3(4,1,-5);
                scene.cylinders[1].direction = normalize(float3(1,4,1));
                scene.cylinders[1].radius = .4;
                scene.cylinders[1].material = getPlasticMaterial();
                scene.cylinders[1].enable = true;
                //getGlassMaterial();

                // compute color for fragment
                //float4 fragCoord = i.vertex;
                //fixed4 col = fixed4(i.uv.x, i.uv.y, 0, 1);

                fixed4 col = fixed4(0, 0, 0, 0);
                //col.r = fragCoord.x;
                //col.g = 0;
                //col.b = 0;
                col.rgb = tonemap(colorForFragment(scene, i.vertex.xy));
                col.a = 1.;
                
                // fixed4 col = tex2D(_MainTex, i.uv);
                return col;
            }
            ENDCG
        }
    }
}
