Shader "Hidden/NewImageEffectShader"
{
    Properties
    {
        _MainTex("Texture 1", 2D) = "white" {}
        _SubTex ("Texture 2", 2D) = "white" {}
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

            #include "UnityCG.cginc"

            struct Ray 
            {
                float3 o;
                float3 d;
            };

            struct Mat {
                Ray r;
            };

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex  : SV_POSITION;
            };

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                //o.screenPosition = ComputeScreenPos(o.vertex);
                o.uv = v.uv;
                return o;
            }

            sampler2D _MainTex;
            sampler2D _SubTex;

            fixed4 frag(v2f i) : SV_Target
            {
                //float2 textureCoordinate = i.screenPosition.xy / i.screenPosition.w;
                Mat c, d;
                Ray a; a.o = float3(1,0,0); a.d = float3(0,1,0);
                Ray b; b.o = float3(0,0,0); b.d = float3(0,0,1);
                c.r = a; d.r = b;
                c = d;
                d.r.o.y = 1.;
                fixed4 col = tex2D(_SubTex, i.uv);
                col.xyz = c.r.o;
                return col;
            }
            ENDCG
        }
    }
}
