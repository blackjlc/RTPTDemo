Shader "Custom/PT"
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

            #define M_PI 3.1415

            static const int sphereCount=4;
            static const int planeCount=4;
            static const int emittingSphereCount=2;
            static const int maxPathLength=3;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };
                
            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            struct Material{
                float lightIntensity;
                float3 lightColor;
                float3 diffuse;
                float3 specular;
                float glossiness;
            };

            struct Sphere{
                float3 position;
                float radius;
                Material material;
            };

            struct Plane{
                float3 normal;
                float d;
                Material material;
            };

            struct Scene{
                Sphere spheres[sphereCount];
                Plane planes[planeCount];
            };

            struct Ray{
                float3 origin;
                float3 direction;
            };
            Ray make_Ray(float3 origin, float3 direction){
                Ray r; r.origin = origin; r.direction = direction; return r;
            }

            // Contains all information pertaining to a ray/object intersection
            struct HitInfo{
            bool hit;
            float t;
            float3 position;
            float3 normal;
            Material material;
            };

            HitInfo make_HitInfo(bool hit, float t, float3 position, float3 normal, Material material){
                HitInfo hf;
                hf.hit = hit;
                hf.t = t;
                hf.position = position;
                hf.normal = normal;
                hf.material = material;
                return hf;
            }

            HitInfo getEmptyHit(){
            Material emptyMaterial;
            emptyMaterial.lightIntensity=0.;
            emptyMaterial.lightColor=float3(0,0,0);
            emptyMaterial.diffuse=float3(0,0,0);
            emptyMaterial.specular=float3(0,0,0);
            emptyMaterial.glossiness=0.;

            HitInfo empHit;
            empHit.hit = false;
            empHit.t = 0;
            empHit.position = float3(0,0,0);
            empHit.normal = float3(0,0,0);
            empHit.material = emptyMaterial;

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
            bool isTInInterval(const float t,const float tMin,const float tMax){
            return t>tMin&&t<tMax;
            }

            // Get the smallest t in an interval
            bool getSmallestTInInterval(float t0,float t1,const float tMin,const float tMax,inout float smallestTInInterval){
            
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
            
            // None was
            return false;
            }

            HitInfo intersectSphere(const Ray ray,const Sphere sphere,const float tMin,const float tMax){
            
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
                
                float3 hitPosition=ray.origin+smallestTInInterval*ray.direction;
                
                float3 normal=
                length(ray.origin-sphere.position)<sphere.radius+.001?
                -normalize(hitPosition-sphere.position):
                normalize(hitPosition-sphere.position);
                
                return make_HitInfo(
                true,
                smallestTInInterval,
                hitPosition,
                normal,
                sphere.material);
            }
            return getEmptyHit();
            }

            HitInfo intersectPlane(Ray ray,Plane plane){
            float t=-(dot(ray.origin,plane.normal)+plane.d)/dot(ray.direction,plane.normal);
            float3 hitPosition=ray.origin+t*ray.direction;
            return make_HitInfo(
                true,
                t,
                hitPosition,
                normalize(plane.normal),
            plane.material);
            return getEmptyHit();
            }

            float lengthSquared(const float3 x){
            return dot(x,x);
            }

            HitInfo intersectScene(const Scene scene,Ray ray,const float tMin,const float tMax)
            {
            HitInfo best_hit_info;
            best_hit_info.t=tMax;
            best_hit_info.hit=false;
            
            for(int i=0;i<sphereCount;++i){
                Sphere sphere=scene.spheres[i];
                HitInfo hit_info=intersectSphere(ray,sphere,tMin,tMax);
                
                if(hit_info.hit&&
                hit_info.t<best_hit_info.t&&
                hit_info.t>tMin)
                {
                best_hit_info=hit_info;
                }
            }
            
            for(int i=0;i<planeCount;++i){
                Plane plane=scene.planes[i];
                HitInfo hit_info=intersectPlane(ray,plane);
                
                if(hit_info.hit&&
                hit_info.t<best_hit_info.t&&
                hit_info.t>tMin)
                {
                best_hit_info=hit_info;
                }
            }
            
            return best_hit_info;
            }

            // Converts a random integer in 15 bits to a float in (0, 1)
            float randomInetegerToRandomFloat(int i){
            return float(i)/32768.;
            }

            // Returns a random integer for every pixel and dimension that remains the same in all iterations
            int pixelIntegerSeed(const int dimensionIndex){
            float3 p=float3(_ScreenParams.xy,dimensionIndex);
            float3 r=float3(23.14069263277926,2.665144142690225,7.358926345);
            return int(32768.*frac(cos(dot(p,r))*123456.));
            }

            // Returns a random float for every pixel that remains the same in all iterations
            float pixelSeed(const int dimensionIndex){
            return randomInetegerToRandomFloat(pixelIntegerSeed(dimensionIndex));
            }

            // The global random seed of this iteration
            // It will be set to a new random value in each step
            uniform int globalSeed;
            int randomSeed;
            void initRandomSequence(){
            randomSeed=globalSeed+pixelIntegerSeed(0);
            }

            // Computesinteger  x modulo y not available in most WEBGL SL implementations
            int mod(const int x,const int y){
            return int(float(x)-floor(float(x)/float(y))*float(y));
            }

            // Returns the next integer in a pseudo-random sequence
            int rand(){
            randomSeed=randomSeed*1103515245+12345;
            return mod(randomSeed/65536,32768);
            }

            // Returns the next float in this pixels pseudo-random sequence
            float uniformRandom(){
            return randomInetegerToRandomFloat(rand());
            }

            // Returns the ith prime number for the first 20
            const int maxDimensionCount=10;
            int prime(const int index){
            if(index==0)return 2;
            if(index==1)return 3;
            if(index==2)return 5;
            if(index==3)return 7;
            if(index==4)return 11;
            if(index==5)return 13;
            if(index==6)return 17;
            if(index==7)return 19;
            if(index==8)return 23;
            if(index==9)return 29;
            if(index==10)return 31;
            if(index==11)return 37;
            if(index==12)return 41;
            if(index==13)return 43;
            if(index==14)return 47;
            if(index==15)return 53;
            return 2;
            }

            float halton(const int sampleIndex,const int dimensionIndex){
                // Find the base with the dimension index smaller than the maximum count
                int d=prime(mod(dimensionIndex,maxDimensionCount));
                float f=1.;
                int j=sampleIndex;
                float result=0.;
                // Keep running as long as j > 0
                for(int i=0;i<100000;i++){
                    if(j>0){
                        // Generate numbers in the sequence
                        f/=float(d);
                        result+=f*float(mod(j,d));
                        j=int(float(j)/float(d));
                    }
                    else{
                        // Use a random float per-pixel and per-dimension to offset the sample
                        result+=pixelSeed(dimensionIndex);
                        // Make sure the sample is inside the pixel
                        if(result>1.)
                            result-=1.;
                        return result;
                    }
                }
                return result;
            }

            // This is the index of the sample controlled by the framework.
            // It increments by one in every call of this shader
            uniform int baseSampleIndex;

            // Returns a well-distributed number in (0,1) for the dimension dimensionIndex
            float sample(const int dimensionIndex){
            return halton(baseSampleIndex,dimensionIndex);
            }

            // This is a helper function to sample two-dimensionaly in dimension dimensionIndex
            float2 sample2(const int dimensionIndex){
            return float2(sample(dimensionIndex+0),sample(dimensionIndex+1));
            }

            float3 sample3(const int dimensionIndex){
            return float3(sample(dimensionIndex+0),sample(dimensionIndex+1),sample(dimensionIndex+2));
            }

            // This is a register of all dimensions that we will want to sample.
            // Thanks to Iliyan Georgiev from Solid Angle for explaining proper housekeeping of sample dimensions in ranomdized Quasi-Monte Carlo
            //
            // So if we want to use lens sampling, we call sample(LENS_SAMPLE_DIMENSION).
            //
            // There are infinitely many path sampling dimensions.
            // These start at PATH_SAMPLE_DIMENSION.
            // The 2D sample pair for vertex i is at PATH_SAMPLE_DIMENSION + PATH_SAMPLE_DIMENSION_MULTIPLIER * i + 0
            #define ANTI_ALIAS_SAMPLE_DIMENSION 0
            #define LENS_SAMPLE_DIMENSION 2
            #define PATH_SAMPLE_DIMENSION 4

            // This is 2 for two dimensions and 2 as we use it for two purposese: NEE and path connection
            #define PATH_SAMPLE_DIMENSION_MULTIPLIER 2

            float3 randomDirection(const int dimensionIndex){
                int sampleDimension=PATH_SAMPLE_DIMENSION+PATH_SAMPLE_DIMENSION_MULTIPLIER*dimensionIndex;
                float2 ksi=sample2(sampleDimension);
                //float ksi0 = sample(sampleDimension);
                //float ksi1 = sample(sampleDimension);
                float theta=acos(2.*ksi.x-1.);
                float phi=ksi.y*2.*M_PI;
                return float3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
            }

            float3 getEmission(const Material material,const float3 normal){
                if(material.lightIntensity>0.)
                    return material.lightIntensity*material.lightColor;
                else return material.diffuse;
            }

            float3 getReflectance(
            const Material material,
            const float3 normal,
            const float3 inDirection,
            const float3 outDirection)
            {
            // Compute the Physically-correct Phong, specular * (glossiness + 2)/(2pi)
            // * (cos(angle between the outgoing light and the perfect specular reflective direction))^glossiness
            // make sure the angle <= pi/2 with max()
            return material.specular*(material.glossiness+2.)/(2.*M_PI)
            *pow(max(0.,dot(normalize(outDirection),normalize(reflect(inDirection,normal)))),material.glossiness)
            +material.diffuse/M_PI;
            }

            float3 getGeometricTerm(
            const Material material,
            const float3 normal,
            const float3 inDirection,
            const float3 outDirection)
            {
            // Find the cos angle between the incoming light direction and the surface normal, which is just the dot product between the two normalized floattors.
            // make sure the angle <= pi/2 with max()
            return max(0.,dot(normalize(inDirection),normalize(normal)))*float3(1,1,1);
            }

            float4x4 rotationMatrixFromAngleAxis(float angle,float3 axis)
            {
            axis=normalize(axis);
            float s=sin(angle);
            float c=cos(angle);
            float oc=1.-c;
            
            return float4x4(oc*axis.x*axis.x+c,
                oc*axis.x*axis.y-axis.z*s,
                oc*axis.z*axis.x+axis.y*s,0.,
                oc*axis.x*axis.y+axis.z*s,
                oc*axis.y*axis.y+c,
                oc*axis.y*axis.z-axis.x*s,0.,
                oc*axis.z*axis.x-axis.y*s,
                oc*axis.y*axis.z+axis.x*s,
                oc*axis.z*axis.z+c,0.,
                0.,
                0.,
                0.,
            1.);
            }

            float3 getEmitterPosition(const float3 position,const Sphere sphere,const int dimensionIndex){
            // This is a simplified version: Just report the sphere center. Will not do well with visibility.
            //return sphere.position;
            
            // This is the wrong simplified version: Take a random surface point.
            //return sphere.position + randomDirection(dimensionIndex) * sphere.radius;
            
            // Well we stick our fingers in
            // The ground, heave and
            // Turn the world around
            
            // This has three main steps:
            //   1) Make a direction
            //   2) Orient it so it points to the sphere
            //   3) Find a point on the sphere along this direction

            // Step 1) Make a random direction in a cone orientedd along th up-pointing z direction
            // .. the opening angle of a sphere in a certain distance
            float apexAngle=asin(sphere.radius/length(position-sphere.position));

            // The rotation around the z axis
            float phi=sample(dimensionIndex+1)*2.*M_PI;

            // z is the cosine of the angle.
            // We need a random cosine of the angle (which is notthe same as the cosine of a random angle!)
            float z=lerp(1.,cos(apexAngle),sample(dimensionIndex+0));
            float3 alignedDirection=float3(sqrt(1.-z*z)*cos(phi),sqrt(1.-z*z)*sin(phi),z);

            // Step 2) Rotate the z axis-aligned dirction to point into the direction of the sphere
            float3 direction=normalize(sphere.position-position);
            float rotationAngle=acos(dot(direction,float3(0,0,1)));
            float3 rotationAxis=cross(direction,float3(0,0,1));
            float4x4 rotationMatrix=rotationMatrixFromAngleAxis(rotationAngle,rotationAxis);
            float3 worldDirection=(mul(rotationMatrix,float4(alignedDirection,0))).xyz;

            // Step 3) Send a ray. it feels this should be easier, but Tobias does not see it.
            Ray emitterRay;
            emitterRay.origin=position;
            emitterRay.direction=worldDirection;
            return intersectSphere(emitterRay,sphere,.01,100000.).position;
            }

            float3 samplePath(const Scene scene,const Ray initialRay){

            // Initial result is black
            float3 result=float3(0,0,0);

            Ray incomingRay=initialRay;
            float3 throughput=float3(1,1,1);
            for(int i=0;i<maxPathLength;i++){
            HitInfo hitInfo=intersectScene(scene,incomingRay,.001,10000.);

            if(!hitInfo.hit)return result;

            result+=throughput*getEmission(hitInfo.material,hitInfo.normal);

            Ray outgoingRay;
            float3 reflDir=randomDirection(i);
            if(dot(reflDir,hitInfo.normal)<0.)
            return result;

            outgoingRay=make_Ray(hitInfo.position,reflDir);

            throughput*=getReflectance(hitInfo.material,hitInfo.normal,outgoingRay.direction,incomingRay.direction)
            *getGeometricTerm(hitInfo.material,hitInfo.normal,outgoingRay.direction,incomingRay.direction);

            // With importance sampling, this value woudl change
            float probability=1.;
            throughput/=probability;

            incomingRay=outgoingRay;
            }
            return result;
            }

            uniform int2 resolution;
            Ray getFragCoordRay(const float2 fragCoord){

            float sensorDistance=1.;
            float3 origin=float3(0,0,sensorDistance);
            float2 sensorMin=float2(-1,-.5);
            float2 sensorMax=float2(1,.5);
            float2 pixelSize=(sensorMax-sensorMin)/float2(resolution);
            float3 direction=normalize(float3(sensorMin+pixelSize*fragCoord,-sensorDistance));

            float apertureSize=0.;
            float focalPlane=100.;
            float3 sensorPosition=origin+focalPlane*direction;
            origin.xy+=apertureSize*(sample2(LENS_SAMPLE_DIMENSION)-float2(.5,.5));
            direction=normalize(sensorPosition-origin);

            return make_Ray(origin,direction);
            }

            float3 colorForFragment(const Scene scene,const float2 fragCoord){
                initRandomSequence();

                // Offset the x and y  coordinate with a random float between 0 ~ 1
                float x=fragCoord.x+uniformRandom();
                float y=fragCoord.y+uniformRandom();
                // Make sure the results are within 0 ~ 1
                if(x>1.)
                    x-=1.;
                if(y>1.)
                    y-=1.;
                float2 sampleCoord=float2(x,y);
                return samplePath(scene,getFragCoordRay(sampleCoord));
            }

            void loadScene1(inout Scene scene){

            scene.spheres[0].position=float3(6,2,-12);
            scene.spheres[0].radius=2.;
            scene.spheres[0].material.lightIntensity=200.;
            scene.spheres[0].material.lightColor=float3(.9,.5,.3);
            scene.spheres[0].material.diffuse=float3(0.,0,0);
            scene.spheres[0].material.specular=float3(0.,0,0);
            scene.spheres[0].material.glossiness=10.;

            scene.spheres[1].position=float3(-7,-2,-11);
            scene.spheres[1].radius=1.;
            scene.spheres[1].material.lightIntensity=200.;
            scene.spheres[1].material.lightColor=float3(.3,.9,.8);
            scene.spheres[1].material.diffuse=float3(0.,0,0);
            scene.spheres[1].material.specular=float3(0.,0,0);
            scene.spheres[1].material.glossiness=10.;

            scene.spheres[2].position=float3(-2,-2,-15);
            scene.spheres[2].radius=3.;
            scene.spheres[2].material.lightIntensity=0.;
            scene.spheres[2].material.lightColor=float3(0,0,0);
            scene.spheres[2].material.diffuse=float3(.2,.3,.8);
            scene.spheres[2].material.specular=float3(.2,.2,.2);
            scene.spheres[2].material.glossiness=4.;

            scene.spheres[3].position=float3(3,-3.5,-11);
            scene.spheres[3].radius=1.;
            scene.spheres[3].material.lightIntensity=0.;
            scene.spheres[3].material.lightColor=float3(0.,0.,0.);
            scene.spheres[3].material.diffuse=float3(0.,0.,0.);
            scene.spheres[3].material.specular=float3(1.,1,1);
            scene.spheres[3].material.glossiness=100.;

            scene.planes[0].normal=float3(0,1,0);
            scene.planes[0].d=4.5;
            scene.planes[0].material.lightIntensity=0.;
            scene.planes[0].material.lightColor=float3(0,0,0);
            scene.planes[0].material.diffuse=float3(.8,.8,.8);
            scene.planes[0].material.specular=float3(0,0,0);
            scene.planes[0].material.glossiness=50.;

            scene.planes[1].normal=float3(0,0,1);
            scene.planes[1].d=18.5;
            scene.planes[1].material.lightIntensity=0.;
            scene.planes[1].material.lightColor=float3(0,0,0);
            scene.planes[1].material.diffuse=float3(.9,.6,.3);
            scene.planes[1].material.specular=float3(.02,.02,.02);
            scene.planes[1].material.glossiness=3000.;

            scene.planes[2].normal=float3(1,-.3,0);
            scene.planes[2].d=10.;
            scene.planes[2].material.lightIntensity=0.;
            scene.planes[2].material.lightColor=float3(0,0,0);
            scene.planes[2].material.diffuse=float3(.2,.2,.2);
            scene.planes[2].material.specular=float3(.1,.1,.1);
            scene.planes[2].material.glossiness=100.;

            scene.planes[3].normal=float3(-1,-.3,0);
            scene.planes[3].d=10.;
            scene.planes[3].material.lightIntensity=0.;
            scene.planes[3].material.lightColor=float3(0,0,0);
            scene.planes[3].material.diffuse=float3(.2,.2,.2);
            scene.planes[3].material.specular=float3(.1,.1,.1);
            scene.planes[3].material.glossiness=1000.;
            }

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
                loadScene1(scene);

                fixed4 col = fixed4(0, 0, 0, 0);
                col.rgb = colorForFragment(scene,i.vertex.xy);
                col.a = 1.;
                
                // fixed4 col = tex2D(_MainTex, i.uv);
                return col;
            }
            ENDCG
        }
    }
}
