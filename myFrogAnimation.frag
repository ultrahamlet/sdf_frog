const float epsilon = 0.01;
const float pi = 3.14159265359;
const float halfpi = 1.57079632679;
const float twopi = 6.28318530718;

mat3 rotateMat(vec3 p, float angle, vec3 axis){
    vec3 a = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float r = 1.0 - c;
    mat3 m = mat3(
        a.x * a.x * r + c,
        a.y * a.x * r + a.z * s,
        a.z * a.x * r - a.y * s,
        a.x * a.y * r - a.z * s,
        a.y * a.y * r + c,
        a.z * a.y * r + a.x * s,
        a.x * a.z * r + a.y * s,
        a.y * a.z * r - a.x * s,
        a.z * a.z * r + c
    );
    return m;
}

#define LIGHT normalize(vec3(1.0, 1.0, 0.0))

float displacement(vec3 p){
    return sin(20.141592*p.x)*sin(0.141592*p.y)*sin(20.131592*p.y);
}

//----------------------------------------------------------------------------------
//  from SDF Editor
//----------------------------------------------------------------------------------



float pCapsule(float r, float h, vec3 p)
{
    p.y -= clamp( p.y, 0.0, h );
    return length( p ) - r;
}



// bound
float pEllipsoid(vec3 r, vec3 p)
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}



vec3 mTranslation(vec3 inv_translation, vec3 p)
{
    return p + inv_translation;
}

vec3 mRotation(mat3 inv_rotation, vec3 p)
{
    return inv_rotation * p;
}

vec3 mMirror(vec3 normal, float dist, vec3 p)
{
    float d = max(0.0, dot(normal, p) - dist);
    return p - 2.0 * d * normal;
}


float oUnion(float d1, float d2)
{
    return min(d1, d2);
}

float oSubtraction(float d1, float d2)
{
    return max(d1, -d2);
}


float oSmoothUnion(float k, float d1, float d2)
{
    float h = clamp(0.5 + 0.5*(d1 - d2) / k, 0.0, 1.0);
    return mix(d1, d2, h) - k*h*(1.0-h);
}

float oSmoothSubtraction(float k, float d1, float d2)
{
    float h = clamp(0.5 - 0.5*(d1 + d2) / k, 0.0, 1.0);
    return mix(d1, -d2, h) + k*h*(1.0-h);
}


//-----------------------------------------------------------------
mat2 rotate(float a)
{
    float s=sin(a),c=cos(a);
    return mat2(c,s,-s,c);
}

mat3 lookat(in vec3 eye, in vec3 target)
{
	vec3 w = normalize(target-eye);
	vec3 u = normalize(cross(w,vec3(0.0,1.0,0.0)));
	vec3 v = normalize(cross(u,w));
    return mat3(u,v,w);
}

float lengthN(in vec2 p, in float n)
{
	//p = pow(abs(p), vec2(n));
    p = pow(abs(p), vec2(n));
	return pow(p.x+p.y, 1.0/n);
}


mat3 RoIn_1 = mat3(1, 0, 0, 0, -0, 1, 0, -1, -0);
vec3 TrIn_2 = vec3(-0, -0.493172, -3.01572);
mat3 RoIn_3 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_4 = vec3(1.99929, 1.97652, 1.28);
vec3 MiNo_5 = vec3(-1, 0, 0);
float MiDi_6 = 0.;
mat3 RoIn_7 = mat3(1, 0, 0, 0, -0, 1, 0, -1, -0);
vec3 TrIn_8 = vec3(-1.67, -2.32, -0);
mat3 RoIn_9 = mat3(-0, 1, 0, -1, -0, 0, 0, 0, 1);
float CaRa_10 = 0.383695;
float CaHe_11 = 1.6905;
vec3 TrIn_12 = vec3(-0, -2.03, -0);
mat3 RoIn_13 = mat3(0.866025, 0.5, -0, -0.5, 0.866025, 0, 0, -0, 1);
float CaRa_14 = 0.383695;
float CaHe_15 = 1.444452;
vec3 TrIn_16 = vec3(-0, -1.83, -0);
mat3 RoIn_17 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
float CaRa_18 = 0.47;
float CaHe_19 = 0.319634;
vec3 TrIn_20 = vec3(-0, -0.73, -0);
mat3 RoIn_21 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
float CaRa_22 = 0.147221;
float CaHe_23 = 1.21;
vec3 TrIn_24 = vec3(0.25, -0.82, -0);
mat3 RoIn_25 = mat3(0.866025, -0.5, 0, 0.5, 0.866025, -0, -0, 0, 1);
float CaRa_26 = 0.15;
float CaHe_27 = 0.93;
vec3 TrIn_28 = vec3(-0.28, -0.59, -0);
mat3 RoIn_29 = mat3(0.939693, 0.34202, -0, -0.34202, 0.939693, 0, 0, -0, 1);
float CaRa_30 = 0.15;
float CaHe_31 = 0.93;
vec3 TrIn_32 = vec3(-0.13, -0.03, -0);
mat3 RoIn_33 = mat3(0.652231, 0.758021, -0, -0.758021, 0.652231, 0, 0, -0, 1);
float CaRa_34 = 0.15;
float CaHe_35 = 0.93;
mat3 RoIn_36 = mat3(1, 0, 0, 0, -0, 1, 0, -1, -0);
vec3 TrIn_37 = vec3(-0.82, -4.86, -0);
mat3 RoIn_38 = mat3(0.939693, 0.34202, -0, -0.34202, 0.939693, 0, 0, -0, 1);
float CaRa_39 = 0.383695;
float CaHe_40 = 1.6905;
vec3 TrIn_41 = vec3(-0, -2.2, -0);
mat3 RoIn_42 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
float CaRa_43 = 0.383695;
float CaHe_44 = 1.444452;
vec3 TrIn_45 = vec3(-0, -1.99, -0);
mat3 RoIn_46 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
float CaRa_47 = 0.477884;
float CaHe_48 = 0.319634;
vec3 TrIn_49 = vec3(-0, -0.53, -0);
mat3 RoIn_50 = mat3(1, 0, -0, -0, 0.766045, 0.642788, 0, -0.642788, 0.766045);
float CaRa_51 = 0.147221;
float CaHe_52 = 1.13;
vec3 TrIn_53 = vec3(-0.09, -0.51, -0);
mat3 RoIn_54 = mat3(0.866025, 0.5, -0, -0.469846, 0.813798, 0.34202, 0.17101, -0.296198, 0.939693);
float CaRa_55 = 0.15;
float CaHe_56 = 1.;
vec3 TrIn_57 = vec3(-0.32, -0.31, -0);
mat3 RoIn_58 = mat3(0.766045, 0.642788, 0, -0.633022, 0.754407, 0.173648, 0.111619, -0.133022, 0.984808);
float CaRa_59 = 0.15;
float CaHe_60 = 0.93;
vec3 TrIn_61 = vec3(-0.19, -0.49, -0);
mat3 RoIn_62 = mat3(0.505883, -0.862602, 0, 0.862602, 0.505883, -0, -0, 0, 1);
float CaRa_63 = 0.15;
float CaHe_64 = 0.93;
vec3 MiNo_65 = vec3(-1, 0, 0);
float MiDi_66 = 0.;
vec3 TrIn_67 = vec3(-1.23772, -0.467108, -4.0436);
mat3 RoIn_68 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_69 = vec3(0.625824, 0.357147, 0.8914);
vec3 TrIn_70 = vec3(-0, 0.876456, -3.01572);
mat3 RoIn_71 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_72 = vec3(1.33164, 1.29976, 0.572888);
vec3 MiNo_73 = vec3(-1, 0, 0);
float MiDi_74 = 0.;
vec3 TrIn_75 = vec3(-1.186, -0.2179, -4.0436);
mat3 RoIn_76 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_77 = vec3(0.4139, 0.357147, 0.66314);
vec3 MiNo_78 = vec3(-1, 0, 0);
float MiDi_79 = 0.;
vec3 TrIn_80 = vec3(-1.23772, -0.467108, -4.0436);
mat3 RoIn_81 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_82 = vec3(0.39685, 0.357147, 0.507432);
vec3 TrIn_83 = vec3(-0, -0.951508, -0.05504);
mat3 RoIn_84 = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
vec3 ElRa_85 = vec3(1.72792, 1.4785, 2.16554);
float SmTr_86 = 0.05;
float SmTr_87 = 0.05;

float sdf(vec3 p0)
{
	float d1;
	float d2;
	float d3;
	float d4;
	float d5;
	float d6;
	float d7;
	float d8;
	float d9;
	float d10;
	float d11;
	float d12;
	float d13;
	float d14;
	float d15;
	float d16;
	float d17;
	float d18;
	float d19;
	float d20;

	{
		vec3 p1 = mRotation(RoIn_1, p0);
		{
			vec3 p2 = mTranslation(TrIn_2, p1);
			{
				vec3 p3 = mRotation(RoIn_3, p2);
				d1 = pEllipsoid(ElRa_4, p3);
				{
					vec3 p4 = mMirror(MiNo_5, MiDi_6, p3);
					{
						
						vec3 p5 = mRotation(RoIn_7, p4);
						{
							vec3 p6 = mTranslation(TrIn_8, p5);
							{
								mat3 mt = rotateMat(p6,halfpi-sin(iTime)*0.5,vec3(0,0,1));
								vec3 p7 = mRotation(mt, p6);
								//vec3 p7 = mRotation(RoIn_9, p6);
								d2 = pCapsule(CaRa_10, CaHe_11, p7);
								{
									vec3 p8 = mTranslation(TrIn_12, p7);
									{
										mat3 mt = rotateMat(p7,halfpi+sin(iTime)*0.4,vec3(0,0,1));
                                        vec3 p9 = mRotation(mt, p8);
										//vec3 p9 = mRotation(RoIn_13, p8);
										d3 = pCapsule(CaRa_14, CaHe_15, p9);
										{
											vec3 p10 = mTranslation(TrIn_16, p9);
											{
												vec3 p11 = mRotation(RoIn_17, p10);
												d4 = pCapsule(CaRa_18, CaHe_19, p11);
												{
													vec3 p12 = mTranslation(TrIn_20, p11);
													{
														vec3 p13 = mRotation(RoIn_21, p12);
														d5 = pCapsule(CaRa_22, CaHe_23, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_24, p11);
													{
														vec3 p13 = mRotation(RoIn_25, p12);
														d6 = pCapsule(CaRa_26, CaHe_27, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_28, p11);
													{
														vec3 p13 = mRotation(RoIn_29, p12);
														d7 = pCapsule(CaRa_30, CaHe_31, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_32, p11);
													{
														vec3 p13 = mRotation(RoIn_33, p12);
														d8 = pCapsule(CaRa_34, CaHe_35, p13);
													}
												}
											}
										}
									}
								}
							}
						}
					}
					{
						vec3 p5 = mRotation(RoIn_36, p4);
						{
							vec3 p6 = mTranslation(TrIn_37, p5);
							{
								mat3 mt = rotateMat(p6,halfpi*0.5+cos(iTime)*0.5,vec3(0,0,1));
								vec3 p7 = mRotation(mt, p6);
								//vec3 p7 = mRotation(RoIn_38, p6);
								d9 = pCapsule(CaRa_39, CaHe_40, p7);
								{
									vec3 p8 = mTranslation(TrIn_41, p7);
									{
										mat3 mt = rotateMat(p8,-halfpi*0.2+sin(iTime)*0.4,vec3(0,0,1));
										vec3 p9 = mRotation(mt, p8);
										//vec3 p9 = mRotation(RoIn_42, p8);
										d10 = pCapsule(CaRa_43, CaHe_44, p9);
										{
											vec3 p10 = mTranslation(TrIn_45, p9);
											{
												vec3 p11 = mRotation(RoIn_46, p10);
												d11 = pCapsule(CaRa_47, CaHe_48, p11);
												{
													vec3 p12 = mTranslation(TrIn_49, p11);
													{
														vec3 p13 = mRotation(RoIn_50, p12);
														d12 = pCapsule(CaRa_51, CaHe_52, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_53, p11);
													{
														vec3 p13 = mRotation(RoIn_54, p12);
														d13 = pCapsule(CaRa_55, CaHe_56, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_57, p11);
													{
														vec3 p13 = mRotation(RoIn_58, p12);
														d14 = pCapsule(CaRa_59, CaHe_60, p13);
													}
												}
												{
													vec3 p12 = mTranslation(TrIn_61, p11);
													{
														vec3 p13 = mRotation(RoIn_62, p12);
														d15 = pCapsule(CaRa_63, CaHe_64, p13);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		{
			vec3 p2 = mMirror(MiNo_65, MiDi_66, p1);
			{
				vec3 p3 = mTranslation(TrIn_67, p2);
				{
					vec3 p4 = mRotation(RoIn_68, p3);
					d16 = pEllipsoid(ElRa_69, p4);
				}
			}
		}
		{
			vec3 p2 = mTranslation(TrIn_70, p1);
			{
				vec3 p3 = mRotation(RoIn_71, p2);
				d17 = pEllipsoid(ElRa_72, p3);
			}
		}
		{
			vec3 p2 = mMirror(MiNo_73, MiDi_74, p1);
			{
				vec3 p3 = mTranslation(TrIn_75, p2);
				{
					vec3 p4 = mRotation(RoIn_76, p3);
					d18 = pEllipsoid(ElRa_77, p4);
				}
			}
		}
		{
			vec3 p2 = mMirror(MiNo_78, MiDi_79, p1);
			{
				vec3 p3 = mTranslation(TrIn_80, p2);
				{
					vec3 p4 = mRotation(RoIn_81, p3);
					d19 = pEllipsoid(ElRa_82, p4);
				}
			}
		}
		{
			vec3 p2 = mTranslation(TrIn_83, p1);
			{
				vec3 p3 = mRotation(RoIn_84, p2);
				d20 = pEllipsoid(ElRa_85, p3);
			}
		}
	}
	return oUnion(oSmoothSubtraction(SmTr_86, d1, d17), oUnion(oSmoothSubtraction(SmTr_87, d16, d18), oUnion(d19, oUnion(d20, oUnion(d5, oUnion(d4, oUnion(d3, oUnion(d2, oUnion(d6, oUnion(d8, oUnion(d7, oUnion(d12, oUnion(d11, oUnion(d10, oUnion(d9, oUnion(d13, oUnion(d14, d15)))))))))))))))));
}
#define ZERO (min(iFrame,0))


vec3 calcNormal( in vec3 pos)
{

    vec2 e = vec2(1.0,-1.0)*0.5773*0.001;
    return normalize( e.xyy*sdf( pos + e.xyy) + 
					  e.yyx*sdf( pos + e.yyx) + 
					  e.yxy*sdf( pos + e.yxy) + 
					  e.xxx*sdf( pos + e.xxx) );

}

float calcOcclusion( in vec3 pos, in vec3 nor)
{
    float res = 0.0;
    vec3 aopos = pos;
    
    for( int i=0; i<4; i++ )
    {   
        aopos = pos + nor*0.2*float(i);
        float d = sdf(aopos);
        res += d;
    }

    return clamp(res, 0.0, 1.0);   
}



float castRay(in vec3 ro, vec3 rd)
{
   float t = 0.0;
    for(int i = 0; i < 100; i++ )
    {
        vec3 pos = ro+t*rd;
        float h = sdf(pos);
       
        if(h < 0.001) break;
            t += h;
        if(t > 20.) break;
        
    }
    if(t > 20.0) t = -1.0;
    return t;

}


vec3 env_color(vec3 dir)
{
	if (dir.y > 0.0)
		return mix(vec3(0.0, 0.5, 1.0), vec3(0.0, 0.1, 0.8), dir.y);
	else
		return mix(vec3(0.0, 0.5, 1.0), vec3(0.8, 0.7, 0.6), pow(-dir.y, 0.5));
}



void mainImage(out vec4 fragColor, in vec2 fragCoord)
{

    vec2 p = (2.0*fragCoord -iResolution.xy)/iResolution.y;
   

    float fi = 10.*iMouse.x/iResolution.x;
    float th = 10.*iMouse.y/iResolution.y;

    vec3 ro = vec3(5.0*cos(fi)*sin(th),5.0*sin(th)*sin(fi),15.0*cos(th));
    
    vec3 ta = vec3(0.0,0.0,0.0);
    vec3 ww = normalize(ta-ro);
    vec3 uu = normalize(cross(ww,vec3(0,1,0)));
    vec3 vv = normalize(cross(uu,ww));
    vec3 rd = normalize(p.x*uu+p.y*vv+1.5*ww);
    //vec3 col = vec3(0.4,0.75,0.9) - 0.5*rd.y;
    vec3 col = vec3(0.25,0.72,0.75) - 0.005*rd.y;
    //col = mix(col, vec3(0.7,0.75,0.8),exp(-10.0*rd.y));
    float t = castRay(ro,rd);
    
    if(t > 0.0){
       
        col = vec3(1.0);
        vec3 pos = ro + t *  rd;
        vec3 nor = calcNormal(pos);
        //float focc = res.w;
        //float focc = 2.5;  //  focc = clamp(4.0*res.z,0.0,1.0); -->candy
        float occ = calcOcclusion(pos, nor);//*focc;
        vec3 mate = vec3(0.18);
        vec3 sun_dir = normalize(vec3(0.8,0.4,0.2));
        float sun_sha = step(castRay(pos+nor*0.001, sun_dir), 0.0);
        float sun_dif = clamp(dot(nor,sun_dir),0.0,1.0);
        float sky_dif = clamp(0.5+0.5*dot(nor,vec3(0.0,1.0,0.0)),0.0,1.0);
        float bou_dif = clamp(0.5+0.5*dot(nor,vec3(0.0,-1.0,0.0)),0.0,1.0);
        col = mate*vec3(7.0,4.5,3.0)*sun_dif*sun_sha;
        col += mate*vec3(0.5,0.8,0.9) *sky_dif*occ;
        col += mate*vec3(0.7,0.3,0.2) *bou_dif*occ;

    } else {
        col = env_color(rd);
        //col = mix(col, vec3(0.7,0.75,0.9),exp(-10.0*rd.y));
    }
    col = pow(col,vec3(0.4545));
    fragColor = vec4(col,1.0);

}

