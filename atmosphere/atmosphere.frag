precision mediump float;

const float PI = 3.14159;
int step = 32;
float cameraFOV = 1.5;
float cameraHeight = 1.0; // km

// Radii of Earth and its atmosphere in kilometers
const float earthRadius = 6371.;
const float atmosphereRadius = 6471.;

// Rayleigh phase function (scattering by molecules in atmosphere)
float phaseRayleigh(float mu) {
	return (3. * (1. + mu*mu)) / (16. * PI);
}

// Henyey-Greenstein phase function (scattering by larger particles like aerosols or water droplets)
float phaseHG(float mu, float g) {
	float numerator = 1. - g*g;
	float denominator = 1. + g*g - 2.*g*mu;
	return numerator / (4. * PI * denominator * sqrt(denominator));
}

// Density of Rayleigh scatterers at a given height
vec3 densityRayleigh(float height) {
	return vec3(3.8e-3, 1.36e-2, 3.35e-2) * exp(-height / 8.0);
}

// Density of Mie scatterers at a given height
vec3 densityMie(float height) {
	return vec3(2.1e-2) * exp(-height / 1.2);
}

// Calculate intersections of a ray with a sphere of given radius
vec2 raySphereIntersect(vec3 rOrigin, vec3 rDir, float radius) {
	float rOrD = dot(rOrigin, rDir);
	float rDirSquared = dot(rDir, rDir);
	float det = rOrD * rOrD - (dot(rOrigin, rOrigin) - radius*radius) * rDirSquared;
	if (det < 0.) {
		return vec2(-1.);
	}
	float sqDet = sqrt(det);
	return (-rOrD + vec2(-sqDet, sqDet)) / rDirSquared;
}

void main() {
	// Convert screen and mouse coordinates to UV coordinates in range of [-1, 1]
	vec2 screenUV = (2.*gl_FragCoord.xy - iResolution.xy) / max(iResolution.x, iResolution.y);
	vec2 mouseUV = (2.*iMouse.xy - iResolution.xy) / max(iResolution.x, iResolution.y);

	// Calculate ray origin and direction
	vec3 rOrigin = vec3(0., 0., earthRadius + cameraHeight);
	vec3 rDir = normalize(vec3(1, screenUV.x * cameraFOV, screenUV.y * cameraFOV));
	vec3 sunDir = normalize(vec3(1, mouseUV.x * cameraFOV, mouseUV.y * cameraFOV));

	// Initialize color accumulator
	vec3 color = vec3(0);

	// Raytracing for atmosphere scattering
	vec2 rAtmosphereRange = raySphereIntersect(rOrigin, rDir, atmosphereRadius);
	vec2 rEarthRange = raySphereIntersect(rOrigin, rDir, earthRadius);

	// If there's no intersection with the atmosphere, the pixel is set to black
	if (rAtmosphereRange.y <= 0.) {
		gl_FragColor = vec4(0, 0, 0, 1);
		return;
	}

	// Determine the range of the ray within the atmosphere
	vec2 rRange;
	if(rEarthRange.x > 0. ) {
		rRange = vec2(max(0., rAtmosphereRange.x), rEarthRange.x);
	}
	else {
		rRange = vec2(max(0., rAtmosphereRange.x), rAtmosphereRange.y);
	}

	// Ray step for traversal through atmosphere
	vec3 rStep = (rRange.y - rRange.x) / float(step) * rDir;
	vec3 rPosition = rOrigin + rRange.x * rDir;
	vec3 rDepth = vec3(0);
	float mu = dot(sunDir, rDir);

	// Traverse the atmosphere with ray
	for (int i = 0; i < step; i++) {
		rPosition += rStep / 2.;

		float rHeight = length(rPosition) - earthRadius;
		vec3 densityScatter = densityMie(rHeight) + densityRayleigh(rHeight);
		vec3 inScatter = densityMie(rHeight) * phaseHG(mu, 0.76) + densityRayleigh(rHeight) * phaseRayleigh(mu);
		rDepth += densityScatter * length(rStep) / 2.;

		// In-scattering ray
		vec2 r1AtmosphereRange = raySphereIntersect(rPosition, sunDir, atmosphereRadius);
		vec2 r1EarthRange = raySphereIntersect(rPosition, sunDir, earthRadius);

		// Continue if ray hits the ground
		if (r1EarthRange.y > 0.) {
			continue;
		}

		// Determine the range of the in-scattering ray within the atmosphere
		vec2 r1Range = vec2(max(0., r1AtmosphereRange.x), r1AtmosphereRange.y);
		vec3 r1Step = (r1Range.y - r1Range.x) / float(step) * sunDir;
		vec3 r1Position = rPosition;
		vec3 r1Depth = vec3(0);

		// Traverse the atmosphere with in-scattering ray
		for (int j = 0; j < step; j++) {
			r1Position += r1Step / 2.;

			float r1Height = length(r1Position) - earthRadius;
			vec3 r1Density = densityMie(r1Height) + densityRayleigh(r1Height);

			r1Depth += r1Density * length(r1Step);
			r1Position += r1Step / 2.;
		}

		// Accumulate color based on in-scattering and extinction
		color += exp(-(rDepth + r1Depth)) * inScatter * length(rStep);

		// Update depth and position for next iteration
		rDepth += densityScatter * length(rStep);
		rPosition += rStep / 2.;
	}

	// Apply gamma correction and set fragment color
	gl_FragColor = vec4(pow(color * 20., vec3(1./ 2.2)), 1.);
}
