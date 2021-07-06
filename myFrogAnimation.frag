// Processed by 'GLSL Shader Shrinker' (Shrunk by 3,177 characters)
// (https://github.com/deanthecoder/GLSLShaderShrinker)

mat3 rotateMat(vec3 p, float angle, vec3 axis) {
	vec3 a = normalize(axis);
	float s = sin(angle),
	      c = cos(angle),
	      r = 1. - c;
	return mat3(a.x * a.x * r + c, a.y * a.x * r + a.z * s, a.z * a.x * r - a.y * s, a.x * a.y * r - a.z * s, a.y * a.y * r + c, a.z * a.y * r + a.x * s, a.x * a.z * r + a.y * s, a.y * a.z * r - a.x * s, a.z * a.z * r + c);
}

float pCapsule(float r, float h, vec3 p) {
	p.y -= clamp(p.y, 0., h);
	return length(p) - r;
}

float pEllipsoid(vec3 r, vec3 p) {
	float k0 = length(p / r);
	return k0 * (k0 - 1.) / length(p / (r * r));
}

vec3 mTranslation(vec3 inv_translation, vec3 p) { return p + inv_translation; }

vec3 mRotation(mat3 inv_rotation, vec3 p) { return inv_rotation * p; }

vec3 mMirror(vec3 normal, float dist, vec3 p) { return p - 2. * max(0., dot(normal, p) - dist) * normal; }

float oUnion(float d1, float d2) { return min(d1, d2); }

float oSmoothSubtraction(float k, float d1, float d2) {
	float h = clamp(.5 - .5 * (d1 + d2) / k, 0., 1.);
	return mix(d1, -d2, h) + k * h * (1. - h);
}

float sdf(vec3 p0) {
	float d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20;
	{
		vec3 p1 = mRotation(mat3(1, 0, 0, 0, -0, 1, 0, -1, -0), p0);
		{
			vec3 p2 = mTranslation(vec3(-0, -.493172, -3.01572), p1);
			{
				vec3 p3 = mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p2);
				d1 = pEllipsoid(vec3(1.99929, 1.97652, 1.28), p3);
				{
					vec3 p4 = mMirror(vec3(-1, 0, 0), 0., p3);
					{
						vec3 p5 = mRotation(mat3(1, 0, 0, 0, -0, 1, 0, -1, -0), p4);
						{
							vec3 p6 = mTranslation(vec3(-1.67, -2.32, -0), p5);
							{
								mat3 mt = rotateMat(p6, 1.5707963 - sin(iTime) * .5, vec3(0, 0, 1));
								vec3 p7 = mRotation(mt, p6);
								d2 = pCapsule(.383695, 1.6905, p7);
								{
									vec3 p8 = mTranslation(vec3(-0, -2.03, -0), p7);
									{
										mat3 mt = rotateMat(p7, 1.5707963 + sin(iTime) * .4, vec3(0, 0, 1));
										vec3 p9 = mRotation(mt, p8);
										d3 = pCapsule(.383695, 1.444452, p9);
										{
											vec3 p10 = mTranslation(vec3(-0, -1.83, -0), p9);
											{
												vec3 p11 = mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p10);
												d4 = pCapsule(.47, .319634, p11);
												{
													vec3 p12 = mTranslation(vec3(-0, -.73, -0), p11);
													{ d5 = pCapsule(.147221, 1.21, mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(.25, -.82, -0), p11);
													{ d6 = pCapsule(.15, .93, mRotation(mat3(.866025, -.5, 0, .5, .866025, -0, -0, 0, 1), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(-.28, -.59, -0), p11);
													{ d7 = pCapsule(.15, .93, mRotation(mat3(.939693, .34202, -0, -.34202, .939693, 0, 0, -0, 1), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(-.13, -.03, -0), p11);
													{ d8 = pCapsule(.15, .93, mRotation(mat3(.652231, .758021, -0, -.758021, .652231, 0, 0, -0, 1), p12)); }
												}
											}
										}
									}
								}
							}
						}
					}
					{
						vec3 p5 = mRotation(mat3(1, 0, 0, 0, -0, 1, 0, -1, -0), p4);
						{
							vec3 p6 = mTranslation(vec3(-.82, -4.86, -0), p5);
							{
								mat3 mt = rotateMat(p6, .7854 + cos(iTime) * .5, vec3(0, 0, 1));
								vec3 p7 = mRotation(mt, p6);
								d9 = pCapsule(.383695, 1.6905, p7);
								{
									vec3 p8 = mTranslation(vec3(-0, -2.2, -0), p7);
									{
										mat3 mt = rotateMat(p8, -.31416 + sin(iTime) * .4, vec3(0, 0, 1));
										vec3 p9 = mRotation(mt, p8);
										d10 = pCapsule(.383695, 1.444452, p9);
										{
											vec3 p10 = mTranslation(vec3(-0, -1.99, -0), p9);
											{
												vec3 p11 = mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p10);
												d11 = pCapsule(.477884, .319634, p11);
												{
													vec3 p12 = mTranslation(vec3(-0, -.53, -0), p11);
													{ d12 = pCapsule(.147221, 1.13, mRotation(mat3(1, 0, -0, -0, .766045, .642788, 0, -.642788, .766045), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(-.09, -.51, -0), p11);
													{ d13 = pCapsule(.15, 1., mRotation(mat3(.866025, .5, -0, -.469846, .813798, .34202, .17101, -.296198, .939693), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(-.32, -.31, -0), p11);
													{ d14 = pCapsule(.15, .93, mRotation(mat3(.766045, .642788, 0, -.633022, .754407, .173648, .111619, -.133022, .984808), p12)); }
												}
												{
													vec3 p12 = mTranslation(vec3(-.19, -.49, -0), p11);
													{ d15 = pCapsule(.15, .93, mRotation(mat3(.505883, -.862602, 0, .862602, .505883, -0, -0, 0, 1), p12)); }
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
			vec3 p2 = mMirror(vec3(-1, 0, 0), 0., p1);
			{
				vec3 p3 = mTranslation(vec3(-1.23772, -.467108, -4.0436), p2);
				{ d16 = pEllipsoid(vec3(.625824, .357147, .8914), mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p3)); }
			}
		}
		{
			vec3 p2 = mTranslation(vec3(-0, .876456, -3.01572), p1);
			{ d17 = pEllipsoid(vec3(1.33164, 1.29976, .572888), mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p2)); }
		}
		{
			vec3 p2 = mMirror(vec3(-1, 0, 0), 0., p1);
			{
				vec3 p3 = mTranslation(vec3(-1.186, -.2179, -4.0436), p2);
				{ d18 = pEllipsoid(vec3(.4139, .357147, .66314), mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p3)); }
			}
		}
		{
			vec3 p2 = mMirror(vec3(-1, 0, 0), 0., p1);
			{
				vec3 p3 = mTranslation(vec3(-1.23772, -.467108, -4.0436), p2);
				{ d19 = pEllipsoid(vec3(.39685, .357147, .507432), mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p3)); }
			}
		}
		{
			vec3 p2 = mTranslation(vec3(-0, -.951508, -.05504), p1);
			{ d20 = pEllipsoid(vec3(1.72792, 1.4785, 2.16554), mRotation(mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), p2)); }
		}
	}
	return oUnion(oSmoothSubtraction(.05, d1, d17), oUnion(oSmoothSubtraction(.05, d16, d18), oUnion(d19, oUnion(d20, oUnion(d5, oUnion(d4, oUnion(d3, oUnion(d2, oUnion(d6, oUnion(d8, oUnion(d7, oUnion(d12, oUnion(d11, oUnion(d10, oUnion(d9, oUnion(d13, oUnion(d14, d15)))))))))))))))));
}

vec3 calcNormal(vec3 pos) {
	const vec2 e = vec2(1, -1) * 58e-5;
	return normalize(e.xyy * sdf(pos + e.xyy) + e.yyx * sdf(pos + e.yyx) + e.yxy * sdf(pos + e.yxy) + e.xxx * sdf(pos + e.xxx));
}

float calcOcclusion(vec3 pos, vec3 nor) {
	float res = 0.;
	vec3 aopos = pos;
	for (int i = 0; i < 4; i++) {
		aopos = pos + nor * .2 * float(i);
		float d = sdf(aopos);
		res += d;
	}

	return clamp(res, 0., 1.);
}

float castRay(vec3 ro, vec3 rd) {
	float t = 0.;
	for (int i = 0; i < 100; i++) {
		vec3 pos = ro + t * rd;
		float h = sdf(pos);
		if (h < .001) break;
		t += h;
		if (t > 20.) break;
	}

	if (t > 20.) t = -1.;
	return t;
}

vec3 env_color(vec3 dir) {
	if (dir.y > 0.) return mix(vec3(0, .5, 1), vec3(0, .1, .8), dir.y);
	return mix(vec3(0, .5, 1), vec3(.8, .7, .6), pow(-dir.y, .5));
}

void mainImage(out vec4 fragColor, vec2 fragCoord) {
	vec2 p = (2. * fragCoord - iResolution.xy) / iResolution.y;
	float t,
	      fi = 10. * iMouse.x / iResolution.x,
	      th = 10. * iMouse.y / iResolution.y;
	vec3 ro = vec3(5. * cos(fi) * sin(th), 5. * sin(th) * sin(fi), 15. * cos(th)),
	     ww = normalize(vec3(0) - ro),
	     uu = normalize(cross(ww, vec3(0, 1, 0))),
	     rd = normalize(p.x * uu + p.y * normalize(cross(uu, ww)) + 1.5 * ww),
	     col = vec3(.25, .72, .75) - .005 * rd.y;
	t = castRay(ro, rd);
	if (t > 0.) {
		col = vec3(1);
		vec3 pos = ro + t * rd,
		     nor = calcNormal(pos);
		float occ = calcOcclusion(pos, nor);
		vec3 sun_dir = normalize(vec3(.8, .4, .2));
		float sun_sha = step(castRay(pos + nor * .001, sun_dir), 0.),
		      sun_dif = clamp(dot(nor, sun_dir), 0., 1.),
		      sky_dif = clamp(.5 + .5 * dot(nor, vec3(0, 1, 0)), 0., 1.),
		      bou_dif = clamp(.5 + .5 * dot(nor, vec3(0, -1, 0)), 0., 1.);
		col = vec3(.18) * vec3(7, 4.5, 3) * sun_dif * sun_sha;
		col += vec3(.18) * vec3(.5, .8, .9) * sky_dif * occ;
		col += vec3(.18) * vec3(.7, .3, .2) * bou_dif * occ;
	}
	else col = env_color(rd);

	fragColor = vec4(pow(col, vec3(.4545)), 1);
}