Shader "Custom/RT_test"
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
            #include "RT_instance.cginc"

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
