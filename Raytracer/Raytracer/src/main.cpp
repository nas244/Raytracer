#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/transform.hpp"
#include "CImg.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <vector>
#include <time.h>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Node
{
	std::vector<tinyobj::index_t>* indices;
	std::vector<int>* index;
	Node* child1;
	Node* child2;
	Node* child3;
	Node* child4;
	Node* child5;
	Node* child6;
	Node* child7;
	Node* child8;
	glm::vec3 c0;
	glm::vec3 c1;
	int level;
	int maxLevel;
};


struct camera
{
	const int width = 300;
	const int height = 300;
	float fov = 45;
	glm::vec3 position = { 0.,0.,0. };
	glm::vec3 lookAt = { 0.,0.,-4. };
	glm::vec3 up = { 0.,1.,0. };
	glm::mat4 worldTo = glm::lookAt(position, lookAt, up);
	glm::mat4 toWorld = glm::inverse(worldTo);
	cimg_library::CImg<unsigned char> renderTarget;
	float lens_radius = .015;
	//float lens_radius = 0;
};

struct directionalLight
{
	glm::vec3 direction = { 0.,1.,1 };
	glm::vec3 invdirection = { -0.,-1,-1 };
	glm::vec3 color = { 1.,1.,1. };
	float intensity = 1;
};

struct scene
{
	int max_depth = 6.;
	int rPerPixel = 4;
	glm::vec3 prevcolor;
	float prevref;
};

struct material
{
	glm::vec3 matcolor;
	float ref = 0.;
	glm::vec3 emissive;
	float ior = 0;
	float gangle = 0;
	float roughness = 0;
	float metalness = 0;
};

struct model
{
	std::vector<tinyobj::shape_t>* shapes = new std::vector<tinyobj::shape_t>();
	std::vector<tinyobj::material_t>* materials = new std::vector<tinyobj::material_t>();
	tinyobj::attrib_t* attrib = new tinyobj::attrib_t();
	glm::mat4 transform;
	glm::mat4 invtransform;
	glm::mat4 invtrastransform;
	struct material* mat;
	int modelnum;
	int modeltype;
	glm::vec3 c0;
	glm::vec3 c1;
	glm::vec3 c2;
	float sr;
};

struct hitInfo
{
	glm::vec3 hitPosition;
	glm::vec3 hitNormal;
	struct materials* hitMaterial;
	float hitDistance = 1000000;
	int model;
	bool didHit;
	float u;
	float v;
	float w;
	int index;
	bool inside;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Node* newNode(struct Node* parent)
{
	Node* n = new Node;
	//n-> = value;
	//n->left = NULL;
	//n->right = NULL;
	return n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float max(float a, float b)
{
	if (a < b)
	{
		return b;
	}
	return a;
}

float min(float a, float b)
{
	if (a > b)
	{
		return b;
	}
	return a;
}

bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) {
		x0 = x1 = -0.5 * b / a;
	}
	else {
		float q = (b > 0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}

	return true;
}

struct hitInfo* TestBoxIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo* hit, glm::vec3 c1, glm::vec3 c2, struct model* mod)
{

	float t_min = (c1.x - ro.x) / ray.x;
	float t_max = (c2.x - ro.x) / ray.x;

	if (t_min > t_max)
	{
		float temp = t_min;
		t_min = t_max;
		t_max = temp;
	}

	float ty_min = (c1.y - ro.y) / ray.y;
	float ty_max = (c2.y - ro.y) / ray.y;

	if (ty_min > ty_max)
	{
		float temp = ty_min;
		ty_min = ty_max;
		ty_max = temp;
	}

	if ((t_min > ty_max) || (ty_min > t_max))
	{
		hit->didHit = false;
		return hit;
	}

	if (ty_min > t_min)
	{
		t_min = ty_min;
	}

	if (ty_max < t_max)
	{
		t_max = ty_max;
	}

	float tz_min = (c1.z - ro.z) / ray.z;
	float tz_max = (c2.z - ro.z) / ray.z;

	if (tz_min > tz_max)
	{
		float temp = tz_min;
		tz_min = tz_max;
		tz_max = temp;
	}

	if ((t_min > tz_max) || (tz_min > t_max))
	{
		hit->didHit = false;
		return hit;
	}

	if (tz_min > t_min)
	{
		t_min = tz_min;
	}

	if (tz_max < t_max)
	{
		t_max = tz_max;
	}

	if (t_min < 0 && t_max < 0)
	{
		hit->didHit = false;
		return hit;
	}

	float distance;
	distance = t_min;
	if (t_min < 0)
	{
		distance = t_max;
	}
	hit->didHit = true;
	if (mod->modeltype == 2)
	{
		if (hit->hitDistance > distance)
		{
			hit->hitDistance = distance;
			hit->hitPosition = ro + ray * hit->hitDistance;
			hit->model = mod->modelnum;
			float espi = .001;
			hit->inside = false;
			if (abs(hit->hitPosition.x - c1.x) < espi)
			{
				hit->hitNormal = glm::vec3(-1., 0., 0.);
			}
			if (abs(hit->hitPosition.x - c2.x) < espi)
			{
				hit->hitNormal = glm::vec3(1., 0., 0.);
			}
			if (abs(hit->hitPosition.y - c1.y) < espi)
			{
				hit->hitNormal = glm::vec3(0., -1., 0.);
			}
			if (abs(hit->hitPosition.y - c2.y) < espi)
			{
				hit->hitNormal = glm::vec3(0., 1., 0.);
			}
			if (abs(hit->hitPosition.z - c1.z) < espi)
			{
				hit->hitNormal = glm::vec3(0., 0., -1.);
			}
			if (abs(hit->hitPosition.z - c2.z) < espi)
			{
				hit->hitNormal = glm::vec3(0., 0., 1.);
			}
		}
		
	}
	else
	{
		hit->hitDistance = distance;
	}
	return hit;
}

struct hitInfo* TestSphIntersection(glm::vec3 ro, glm::vec3 ray, struct model* mod, struct hitInfo* hit)
{
	float distance;

	glm::vec3 so = mod->c0;
	float sr = mod->sr;
	//sphere intersect
	float t0, t1; // solutions for t if the ray intersects 
		// analytic solution
	glm::vec3 L = ro - so;
	float a = glm::dot(ray,ray);
	float b = 2 * glm::dot(ray, L);
	float c = glm::dot(L,L) - sr*sr;
	if (!solveQuadratic(a, b, c, t0, t1))
	{
		hit->didHit = false;
		return hit;
	}

	if (t0 > t1) std::swap(t0, t1);

	bool inside = false;

	if (t0 < 0) {
		t0 = t1; // if t0 is negative, let's use t1 instead 
		inside = true;
		if (t0 < 0)
		{
			hit->didHit = false;
			return hit;// both t0 and t1 are negative
		}
	}

	distance = t0;

	if (distance < hit->hitDistance)
	{
		hit->hitDistance = distance;
		hit->hitPosition = ro + ray * hit->hitDistance;
		hit->model = mod->modelnum;
		hit->hitNormal = hit->hitPosition - so  ;

		hit->inside = inside;
		hit->hitNormal = glm::normalize(hit->hitNormal);
		if (inside)
		{
			hit->hitNormal = hit->hitNormal * -1.0f;
		}
		hit->didHit = true;
	}
	

	return hit;
}


// Ray triangle intersection using the Moller-Trumbore algorithm (scratchapixel.com)
struct hitInfo* TestTriIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo* hit, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2 )
{
	glm::vec3 v0v1 = v1 - v0;
	glm::vec3 v0v2 = v2 - v0;
	glm::vec3 pvec = glm::cross(ray, v0v2);
	float det = glm::dot(v0v1, pvec);

	if (fabs(det) < 0.0001)
	{
		hit->didHit = false;
		return hit;
	}

	float invDet = 1.0f / det;

	glm::vec3 tvec = ro - v0;
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0 || u > 1)
	{
		hit->didHit = false;
		return hit;
	}

	glm::vec3 qvec = glm::cross(tvec, v0v1);
	float v = glm::dot(ray, qvec) * invDet;

	if (v < 0 || u + v > 1)
	{
		hit->didHit = false;
		return hit;
	}

	hit->hitDistance = glm::dot(v0v2, qvec) * invDet;

	if (hit->hitDistance < 0.0001)
	{
		hit->didHit = false;
		return hit;
	}
		
	//hitInfo.hitPosition = rayOrigin + rayDirection * hitInfo.hitDistance;
	//hitInfo.hitMaterial = m_material;

	//h.hitNormal = glm::normalize((1.0f - u - v) * n0 + u * n1 + v * n2);

	hit->didHit = true;
	hit->u = u;
	hit->v = v;
	hit->inside = false;
	if (det < 0.0001)
	{
		hit->inside = true;
	}

	//hit.hitNormal = glm::vec3{ 0.,0.,-1. };
	//glm::vec4 tempNormal = m_transform * glm::vec4(interpolatedNormal, 0);
	//hitInfo.hitNormal = { tempNormal.x, tempNormal.y, tempNormal.z };

	return hit;
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Fresnel's equation from scratchapixel
float fresnel(glm::vec3 ray, glm::vec3 normal, float ior)
{
	//ray = ray * -1.f;
	float kr;
	float cosi = glm::dot(ray, normal);
	float etai = 1, etat = ior;
	if (cosi > 0) { std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		kr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1.0f - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		kr = (Rs * Rs + Rp * Rp) / 2.0f;
	}
	return kr;
}

glm::vec3 fresnelSchlick(float cosTheta, glm::vec3 F0)
{
	float mult = pow(1.0 - cosTheta, 5.0);
	return F0 + (glm::vec3(1.0,1.0,1.0) - F0) * mult;
}

float DistributionGGX(glm::vec3 N, glm::vec3 H, float roughness)
{
	float a = roughness * roughness;
	float a2 = a * a;
	float NdotH = max(0.,glm::dot(N, H));
	float NdotH2 = NdotH * NdotH;

	float num = a2;
	float denom = (NdotH2 * (a2 - 1.0) + 1.0);
	denom = 3.14159 * denom * denom;

	return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
	float r = (roughness + 1.0);
	float k = (r * r) / 8.0;

	float num = NdotV;
	float denom = NdotV * (1.0 - k) + k;

	return num / denom;
}

float GeometrySmith(glm::vec3 N, glm::vec3 V, glm::vec3 L, float roughness)
{
	float NdotV = max(0.,glm::dot(N, V));
	float NdotL = max(0.,glm::dot(N, L));
	float ggx2 = GeometrySchlickGGX(NdotV, roughness);
	float ggx1 = GeometrySchlickGGX(NdotL, roughness);

	return ggx1 * ggx2;
}

//Cook Torrence
glm::vec3 CookTorrence(glm::vec3 camDir, glm::vec3 lightDir, glm::vec3 normal, glm::vec3 ray, float roughness, glm::vec3 albedo, float metallic)
{
	glm::vec3 L = glm::normalize(lightDir);
	glm::vec3 N = normal;
	glm::vec3 V = camDir;
	glm::vec3 H = normalize(V + L);


	//glm::vec3 halfAngle = normalize(camDir + lightDir);
	//float NdotV = max(0., glm::dot(normal, camDir));
	float NdotL = max(0., glm::dot(N, L));
	//float NdotH = max(0., glm::dot(normal, L));
	//float VdotH = max(0., dot(camDir, H));

	float standard = 2.5;
	//float F = fresnel(ray, normal, standard);
	/*
	//Fresnel-Schlick
	glm::vec3 F0 = glm::vec3(0.04);
	F0 = glm::mix(F0, albedo, metallic);
	glm::vec3 F = fresnelSchlick(max(0,glm::dot(camDir, H)), F0);

	//Normal Distribution
	float NH2 = pow(NdotH, 2.0);
	float roughness2 = pow(roughness, 2.0);

	float denom = NH2 * roughness2 + (1.0 - NH2);
	float D = roughness2 / (3.14159 * pow(denom, 2.0));

	//Geometric Attentuation
	float g1 = (NdotL * 2.0) / (NdotL + sqrt(roughness2 + (1.0 - roughness2) * pow(NdotL, 2.0)));
	float g2 = (NdotV * 2.0) / (NdotV + sqrt(roughness2 + (1.0 - roughness2) * pow(NdotV, 2.0)));
	float G = g1 * g2;
	*/

	glm::vec3 F0 = glm::vec3(0.04);
	F0 = glm::mix(F0, albedo, metallic);

	float NDF = DistributionGGX(N, H, roughness);

	float G = GeometrySmith(N, V, L, roughness);
	glm::vec3 F = fresnelSchlick(max(dot(H, V), 0.0), F0);

	//Lambertian Diffuse Coefficent
	glm::vec3 kS = F;
	glm::vec3 kD = glm::vec3(1.0) - kS;
	kD *= 1.0 - metallic;

	//Cook-Torrence
	//float denom2 = (3.14159 * NdotV);
	float denominator = 4.0 * max(0.0, glm::dot(N, V)) * max(0., glm::dot(N, L));
	glm::vec3 specular = (G * F * NDF) / denominator;

	//Lambertian Diffuse + Cook-Torrence
	glm::vec3 final = (((kD * albedo) / 3.14159f) + specular) * glm::vec3(1., 1., 1.) * NdotL;

	final = glm::clamp(final, { 0.,0.,0. }, { 1.,1.,1. });

	//if (final.x > 1.)
	//{
		//std::cout << "kill me" << std::endl;
		//Sleep(5000);
	//}

	//ambient
	//glm::vec3 ambient = glm::vec3(0.03) * albedo * .5f;
	//final += ambient;

	return final;
}

bool overlap(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 max, glm::vec3 min)
{
	glm::vec3 boxvertices[8];
	boxvertices[0] = max;
	boxvertices[1] = min;
	boxvertices[2] = glm::vec3{ max.x,max.y,min.z };
	boxvertices[3] = glm::vec3{ max.x,min.y,max.z };
	boxvertices[4] = glm::vec3{ max.x,min.y,min.z };
	boxvertices[5] = glm::vec3{ min.x,max.y,max.z };
	boxvertices[6] = glm::vec3{ min.x,max.y,min.z };
	boxvertices[7] = glm::vec3{ min.x,min.y,max.z };

	glm::vec3 axes[7];
	axes[0] = glm::vec3{ 1.,0.,0. };
	axes[1] = glm::vec3{ 0.,1.,0. };
	axes[2] = glm::vec3{ 0.,0.,1. };

	glm::vec3 v0v1 = v1 - v0;
	glm::vec3 v0v2 = v2 - v0;
	glm::vec3 v1v2 = v2 - v1;

	glm::vec3 normal = glm::cross(v0v1, v0v2);
	normal = glm::normalize(normal);
	
	glm::vec3 axis0 = glm::normalize(glm::cross(v0v1, normal));
	glm::vec3 axis1 = glm::normalize(glm::cross(v0v2, normal));
	glm::vec3 axis2 = glm::normalize(glm::cross(v1v2, normal));

	axes[3] = normal;
	axes[4] = axis0;
	axes[5] = axis1;
	axes[6] = axis2;

	for (int i = 0; i < 7; i++)
	{
		int boundmin = 100000;
		int boundmax = -100000;

		for (int j = 0; j < 8; j++)
		{
			float temp = glm::dot(boxvertices[j], axes[i]);
			if (boundmin > temp)
			{
				boundmin = temp;
			}
			if (boundmax < temp)
			{
				boundmax = temp;
			}
		}

		float trianglemin = 100000;
		float trianglemax = -100000;

		float temp = glm::dot(v0, axes[i]);
		if (trianglemin > temp)
		{
			trianglemin = temp;
		}
		if (trianglemax < temp)
		{
			trianglemax = temp;
		}
		temp = glm::dot(v1, axes[i]);
		if (trianglemin > temp)
		{
			trianglemin = temp;
		}
		if (trianglemax < temp)
		{
			trianglemax = temp;
		}
		temp = glm::dot(v2, axes[i]);
		if (trianglemin > temp)
		{
			trianglemin = temp;
		}
		if (trianglemax < temp)
		{
			trianglemax = temp;
		}

		if (trianglemin > boundmax || trianglemax < boundmin)
		{
			return false;
		}
	}

	return true;
}

struct hitInfo* traverse(struct model* mod, struct Node n, struct hitInfo* hit, glm::vec3 ro, glm::vec3 ray)
{
	struct hitInfo* box = new struct hitInfo;
	box = TestBoxIntersection(ro, ray, box, n.c0, n.c1, mod);

	if (!box->didHit)
	{
		delete box;
		box = NULL;
		return hit;
	}

	

	
	if (n.level == n.maxLevel )
	{
		struct hitInfo* triangle = new struct hitInfo;
		glm::vec3 v0, v1, v2;
		int closest = -1;
		float closestHit = 1000000;
		for (int j = 0; j < n.index->size(); j += 1)
		{
			int index_0 = (*n.indices)[(*n.index)[j]].vertex_index * 3;
			v0.x = mod->attrib->vertices[index_0];
			v0.y = mod->attrib->vertices[index_0 + 1];
			v0.z = mod->attrib->vertices[index_0 + 2];
			int normalIndex_0 = (*n.indices)[(*n.index)[j]].normal_index * 3;
			glm::vec3 normal_0 = { mod->attrib->normals[normalIndex_0],mod->attrib->normals[normalIndex_0 + 1],mod->attrib->normals[normalIndex_0 + 2] };

			int index_1 = (*n.indices)[(*n.index)[j] + 1].vertex_index * 3;
			v1.x = mod->attrib->vertices[index_1];
			v1.y = mod->attrib->vertices[index_1 + 1];
			v1.z = mod->attrib->vertices[index_1 + 2];
			int normalIndex_1 = (*n.indices)[(*n.index)[j] + 1].normal_index * 3;
			glm::vec3 normal_1 = { mod->attrib->normals[normalIndex_1],mod->attrib->normals[normalIndex_1 + 1],mod->attrib->normals[normalIndex_1 + 2] };

			int index_2 = (*n.indices)[(*n.index)[j] + 2].vertex_index * 3;
			v2.x = mod->attrib->vertices[index_2];
			v2.y = mod->attrib->vertices[index_2 + 1];
			v2.z = mod->attrib->vertices[index_2 + 2];
			int normalIndex_2 = (*n.indices)[(*n.index)[j] + 2].normal_index * 3;
			glm::vec3 normal_2 = { mod->attrib->normals[normalIndex_2],mod->attrib->normals[normalIndex_2 + 1],mod->attrib->normals[normalIndex_2 + 2] };

			triangle = TestTriIntersection(ro, ray, triangle, v0, v1, v2);

			if (triangle->didHit)
			{
				if (triangle->hitDistance < hit->hitDistance)
				{
					hit->didHit = true;
					hit->index = j;
					hit->hitDistance = triangle->hitDistance;
					hit->u = triangle->u;
					hit->v = triangle->v;
					hit->w = 1 - (hit->u + hit->v);
					hit->hitNormal = hit->w * normal_0 + hit->u * normal_1 + hit->v * normal_2;
					if (hit->inside)
					{
						hit->hitNormal * -1.f;
					}
					glm::vec4 temp = mod->invtrastransform * glm::vec4{ hit->hitNormal, 0. };
					hit->hitNormal = glm::normalize(glm::vec3{ temp.x,temp.y,temp.z });
					hit->model = mod->modelnum;
					
				}
			}
		}
		delete triangle;
		triangle = NULL;
		delete box;
		box = NULL;
		return hit;
	}

	if (n.child1 != NULL)
	{
		hit = traverse(mod, (*n.child1), hit, ro, ray);
	}

	if (n.child2 != NULL)
	{
		hit = traverse(mod, (*n.child2), hit, ro, ray);
	}
	if (n.child3 != NULL)
	{
		hit = traverse(mod, (*n.child3), hit, ro, ray);
	}
	if (n.child4 != NULL)
	{
		hit = traverse(mod, (*n.child4), hit, ro, ray);
	}
	if (n.child5 != NULL)
	{
		hit = traverse(mod, (*n.child5), hit, ro, ray);
	}
	if (n.child6 != NULL)
	{
		hit = traverse(mod, (*n.child6), hit, ro, ray);
	}
	if (n.child7 != NULL)
	{
		hit = traverse(mod, (*n.child7), hit, ro, ray);
	}
	if (n.child8 != NULL)
	{
		hit = traverse(mod, (*n.child8), hit, ro, ray);
	}

	
	delete box;
	box = NULL;
	return hit;
}

struct Node* subdivide(struct model mod, struct Node n, glm::vec3 c0, glm::vec3 c1)
{
	struct Node* n1 = new struct Node();
	n1->c0 = c0;
	n1->c1 = c1;
	n1->indices = n.indices;
	n1->maxLevel = n.maxLevel;
	n1->level = n.level + 1;

	if (n1->level > n.maxLevel)
	{
		return NULL;
	}

	float max_x, min_x, max_y, min_y, max_z, min_z;
	if (n1->c0.x > n1->c1.x)
	{
		max_x = n1->c0.x;
		min_x = n1->c1.x;
	}
	else
	{
		max_x = n1->c1.x;
		min_x = n1->c0.x;
	}
	if (n1->c0.y > n1->c1.y)
	{
		max_y = n1->c0.y;
		min_y = n1->c1.y;
	}
	else
	{
		max_y = n1->c1.y;
		min_y = n1->c0.y;
	}
	if (n1->c0.z > n1->c1.z)
	{
		max_z = n1->c0.z;
		min_z = n1->c1.z;
	}
	else
	{
		max_z = n1->c1.z;
		min_z = n1->c0.z;
	}

	glm::vec3 v0 = { 0.,0.,0. };
	glm::vec3 v1 = { 0.,0.,0. };
	glm::vec3 v2 = { 0.,0.,0. };
	std::vector<tinyobj::index_t>* indices;
	std::vector<int>* index = new std::vector<int>();
	bool withinbounds = false;

	for (int j = 0; j < n.index->size(); j += 1)
	{
		withinbounds = false;
		int index_0 = (*n.indices)[(*n.index)[j]].vertex_index * 3;
		v0.x = mod.attrib->vertices[index_0];
		v0.y = mod.attrib->vertices[index_0 + 1];
		v0.z = mod.attrib->vertices[index_0 + 2];

		int index_1 = (*n.indices)[(*n.index)[j] + 1].vertex_index * 3;
		v1.x = mod.attrib->vertices[index_1];
		v1.y = mod.attrib->vertices[index_1 + 1];
		v1.z = mod.attrib->vertices[index_1 + 2];

		int index_2 = (*n.indices)[(*n.index)[j] + 2].vertex_index * 3;
		v2.x = mod.attrib->vertices[index_2];
		v2.y = mod.attrib->vertices[index_2 + 1];
		v2.z = mod.attrib->vertices[index_2 + 2];

		//if ((v0.x <= max_x && v0.x >= min_x) && (v0.y <= max_y && v0.y >= min_y) && (v0.z <= max_z && v0.z >= min_z))
		//{
			//withinbounds = true;
		//}
		//if ((v1.x <= max_x && v1.x >= min_x) && (v1.y <= max_y && v1.y >= min_y) && (v1.z <= max_z && v1.z >= min_z))
		//{
			//withinbounds = true;
		//}
		//if ((v2.x <= max_x && v2.x >= min_x) && (v2.y <= max_y && v2.y >= min_y) && (v2.z <= max_z && v2.z >= min_z))
		//{
			//withinbounds = true;
		//}

		withinbounds = overlap(v0, v1, v2, glm::vec3{ max_x,max_y,max_z }, glm::vec3{ min_x,min_y,min_z });
		
		if (withinbounds)
		{
			index->push_back(j * 3);
		}
	}
	
	n1->index = index;

	float average_x = (max_x + min_x) / 2;
	float average_y = (max_y + min_y) / 2;
	float average_z = (max_z + min_z) / 2;
	

	n1->child1 = subdivide(mod, (*n1), glm::vec3{ max_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child2 = subdivide(mod, (*n1), glm::vec3{ max_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child3 = subdivide(mod, (*n1), glm::vec3{ max_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child4 = subdivide(mod, (*n1), glm::vec3{ max_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child5 = subdivide(mod, (*n1), glm::vec3{ min_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child6 = subdivide(mod, (*n1), glm::vec3{ min_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child7 = subdivide(mod, (*n1), glm::vec3{ min_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child8 = subdivide(mod, (*n1), glm::vec3{ min_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });

	return n1;
}

struct Node BVH(struct model mod, struct Node n)
{
	glm::vec3 v0 = { 0.,0.,0. };
	glm::vec3 v1 = { 0.,0.,0. };
	glm::vec3 v2 = { 0.,0.,0. };
	std::vector<tinyobj::index_t>* indices;
	std::vector<int>* index = new std::vector<int>;

	float max_x = -10000;
	float max_y = -10000;
	float max_z = -10000;
	float min_x = 10000;
	float min_y = 10000;
	float min_z = 10000;

	for (int i = 0; i < mod.shapes->size(); i++)
	{
		std::vector<tinyobj::shape_t>* temp = mod.shapes;
		indices = &(*temp)[i].mesh.indices;
		for (int j = 0; j < indices->size(); j += 3)
		{
			int index_0 = (*indices)[j].vertex_index * 3;
			v0.x = mod.attrib->vertices[index_0];
			v0.y = mod.attrib->vertices[index_0 + 1];
			v0.z = mod.attrib->vertices[index_0 + 2];
			if (v0.x > max_x)
			{
				max_x = v0.x;
			}
			if (v0.x < min_x)
			{
				min_x = v0.x;
			}
			if (v0.y > max_y)
			{
				max_y = v0.y;
			}
			if (v0.y < min_y)
			{
				min_y = v0.y;
			}
			if (v0.z > max_z)
			{
				max_z = v0.z;
			}
			if (v0.z < min_z)
			{
				min_z = v0.z;
			}

			int index_1 = (*indices)[j + 1].vertex_index * 3;
			v1.x = mod.attrib->vertices[index_1];
			v1.y = mod.attrib->vertices[index_1 + 1];
			v1.z = mod.attrib->vertices[index_1 + 2];

			if (v1.x > max_x)
			{
				max_x = v1.x;
			}
			if (v1.x < min_x)
			{
				min_x = v1.x;
			}
			if (v1.y > max_y)
			{
				max_y = v1.y;
			}
			if (v1.y < min_y)
			{
				min_y = v1.y;
			}
			if (v1.z > max_z)
			{
				max_z = v1.z;
			}
			if (v1.z < min_z)
			{
				min_z = v1.z;
			}

			int index_2 = (*indices)[j + 2].vertex_index * 3;
			v2.x = mod.attrib->vertices[index_2];
			v2.y = mod.attrib->vertices[index_2 + 1];
			v2.z = mod.attrib->vertices[index_2 + 2];

			if (v2.x > max_x)
			{
				max_x = v2.x;
			}
			if (v2.x < min_x)
			{
				min_x = v2.x;
			}
			if (v2.y > max_y)
			{
				max_y = v2.y;
			}
			if (v2.y < min_y)
			{
				min_y = v2.y;
			}
			if (v2.z > max_z)
			{
				max_z = v2.z;
			}
			if (v2.z < min_z)
			{
				min_z = v2.z;
			}
			index->push_back(j);
		}
	}


	n.c0 = glm::vec3{ min_x,min_y,min_z };
	n.c1 = glm::vec3{ max_x,max_y,max_z };
	n.indices = indices;
	n.maxLevel = 1;
	n.level = 0;
	n.index = index;

	float average_x = (max_x + min_x) / 2;
	float average_y = (max_y + min_y) / 2;
	float average_z = (max_z + min_z) / 2;

	n.child1 = subdivide(mod, n, glm::vec3{ max_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n.child2 = subdivide(mod, n, glm::vec3{ max_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n.child3 = subdivide(mod, n, glm::vec3{ max_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n.child4 = subdivide(mod, n, glm::vec3{ max_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n.child5 = subdivide(mod, n, glm::vec3{ min_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n.child6 = subdivide(mod, n, glm::vec3{ min_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	n.child7 = subdivide(mod, n, glm::vec3{ min_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n.child8 = subdivide(mod, n, glm::vec3{ min_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	

	return n;
}

glm::vec3 Trace(glm::vec3 ray, glm::vec3 ro, int depth, struct camera c, struct scene s, glm::vec3 color, struct directionalLight light, std::vector<struct model> mod, std::vector<struct Node> n)
{
	
	struct hitInfo* hit = new struct hitInfo;
	if (depth > s.max_depth)
	{
		delete hit;
		hit = NULL;
		return color;
	}



	for (int i = 0; i < mod.size(); i++)
	{
		struct model* m;
		
		m = &mod[i];
		
		if (m->modeltype == 0)
		{
			struct Node node;
			node = n[i];
			glm::vec4 temp = m->invtransform * glm::vec4{ ro, 1. };
			glm::vec3 model_ro = glm::vec3{ temp.x,temp.y,temp.z };
			temp = m->invtransform * glm::vec4{ ray, 0. };
			glm::vec3 model_ray = glm::vec3{ temp.x,temp.y,temp.z };

			hit = traverse(m, node, hit, model_ro, model_ray);
		}
		else if (m->modeltype == 1)
		{
			hit = TestSphIntersection(ro, ray, m, hit);
		}
		else if (m->modeltype == 2)
		{
			hit = TestBoxIntersection(ro, ray, hit, m->c0, m->c1, m);
		}
	}
	

	
	
	if (hit->hitDistance == 1000000)
	{
		delete hit;
		hit = NULL;
		return glm::vec3{ 1.,1.,1. };
	}

	
	hit->hitPosition = ro + ray * hit->hitDistance;

	struct hitInfo* occulsion = new struct hitInfo;
	occulsion->didHit = false;
	light.direction = glm::normalize(light.direction);
	
	for (int i = 0; i < mod.size(); i++)
	{
		struct model* x;

		x =&mod[i];

		if (x->modeltype == 0 && x->mat->ior == 0.)
		{
			struct Node node;
			node = n[i];
			glm::vec4 temp = x->invtransform * glm::vec4{ hit->hitPosition + 0.001f * light.direction, 1. };
			glm::vec3 model_ro = glm::vec3{ temp.x,temp.y,temp.z };
			temp = x->invtransform * glm::vec4{ light.direction, 0. };
			glm::vec3 model_ray = glm::vec3{ temp.x,temp.y,temp.z };

			occulsion = traverse(x, node, occulsion, model_ro, model_ray);
			if (occulsion->didHit)
			{
				break;
			}
		}
		else if (x->modeltype == 1 && x->mat->ior == 0.)
		{
			occulsion = TestSphIntersection(hit->hitPosition + 0.001f * light.direction, light.direction, x, occulsion);
			if (occulsion->didHit)
			{
				break;
			}
		}
		else if (x->modeltype == 2 && x->mat->ior == 0.)
		{
			occulsion = TestBoxIntersection(hit->hitPosition + 0.001f * light.direction, light.direction, occulsion, x->c0, x->c1, x);
			if (occulsion->didHit)
			{
				break;
			}
		}
	}

	if (occulsion->didHit)
	{
		delete hit;
		hit = NULL;
		delete occulsion;
		occulsion = NULL;
		return color;
	}
	
	int temp = hit->model;
	struct material mat = *mod[temp].mat;
	
	float diff;
	diff = glm::dot(hit->hitNormal, light.direction);

	if (diff < 0) 
	{
		diff = 0;
	}

	glm::vec3 diffuse = diff * light.intensity * light.color;
	
	glm::vec3 viewdir = normalize(c.position - hit->hitPosition);
	color = CookTorrence(viewdir, light.direction, hit->hitNormal, ray, mat.roughness, mat.matcolor, mat.metalness);

	if (mat.ref > 0)
	{
		s.prevcolor = color;
		s.prevref = mat.ref;
		glm::vec3 reflect_o;
		glm::vec3 reflect = glm::normalize(glm::reflect(ray, hit->hitNormal));
		if (mod[temp].mat->gangle > 0.00001)
		{
			glm::vec3 trans_point = hit->hitPosition + reflect;
			float offset = tan(mod[temp].mat->gangle / 4.0f);
			glm::vec3 new_point = glm::vec3{ ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset };
			trans_point = trans_point + new_point;
			reflect = trans_point - hit->hitPosition;
		}
		reflect_o = hit->hitPosition + 0.001f * reflect;

		glm::vec3 tempcolor = Trace(reflect, reflect_o, depth + 1, c, s, { 0.,0.,0. }, light, mod, n);
		//if (tempcolor.y > 0)
		//{
			//std::cout << "STOP" << std::endl;
			//Sleep(5000);
		//}
		tempcolor = color * (1 - mat.ref) + tempcolor * mat.ref * mat.matcolor;
		color = tempcolor;
		
	}
	else if (mat.ior > 0)
	{
		float kr;
		s.prevref = mat.ref;
		glm::vec3 reflect_o;
		glm::vec3 refract_o;
		glm::vec3 refract_tempcolor = glm::vec3{ 0,0,0 };

		if (hit->inside)
		{
			kr = fresnel(ray, hit->hitNormal * -1.0f, mat.ior);
		}
		else
		{
			kr = fresnel(ray, hit->hitNormal, mat.ior);
		}

		glm::vec3 reflect = glm::normalize(glm::reflect(ray, hit->hitNormal));
		if (mod[temp].mat->gangle > 0.00001)
		{
			glm::vec3 trans_point = hit->hitPosition + reflect;
			float offset = tan(mod[temp].mat->gangle / 4.0f);
			glm::vec3 new_point = glm::vec3{ ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset };
			trans_point = trans_point + new_point;
			reflect = trans_point - hit->hitPosition;
		}
		reflect_o = hit->hitPosition + 0.001f * reflect;
		glm::vec3 reflect_tempcolor = Trace(reflect, reflect_o, depth + 1, c, s, color, light, mod, n);

		if (kr < 1)
		{
			float ior = 1.f / mat.ior;
			glm::vec3 refract;
			if (hit->inside)
			{
				refract = glm::normalize(glm::refract(ray, hit->hitNormal, mat.ior));
			}
			else
			{
				refract = glm::normalize(glm::refract(ray, hit->hitNormal, ior));
			}
			if (mod[temp].mat->gangle > 0.00001)
			{
				glm::vec3 trans_point = hit->hitPosition + refract;
				float offset = tan(mod[temp].mat->gangle / 4.0f);
				glm::vec3 new_point = glm::vec3{ ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset, ((static_cast <float> (rand()) / (static_cast <float> (RAND_MAX))) * 2 - 1) * offset };
				trans_point = trans_point + new_point;
				refract = trans_point - hit->hitPosition;
			}
			refract_o = hit->hitPosition + 0.001f * refract;
			refract_tempcolor = Trace(refract, refract_o, depth + 1, c, s, color, light, mod, n);
		}

		color = reflect_tempcolor * kr + refract_tempcolor * (1 - kr);
	}

	delete occulsion;
	occulsion = NULL;
	delete hit;
	hit = NULL;

	return color;
}




int main(int argc, char* argv[])
{
	struct camera c;
	struct directionalLight light;
	struct scene s;
	struct model mod1;
	struct model mod2;
	struct model mod3;
	struct model mod4;
	struct model mod5;
	struct model mod6;
	struct material* mat1 = new struct material;
	struct material* mat2 = new struct material;
	struct material* mat3 = new struct material;
	struct material* mat4 = new struct material;
	struct material* mat5 = new struct material;
	struct material* mat6 = new struct material;

	bool niaveRien = false;
	std::vector<std::vector<glm::vec3>> buffer;
	buffer.resize(c.width);
	for (int i = 0; i < c.width; ++i)
	{
		buffer[i].resize(c.height, glm::vec3(0.));
	}

	double inversewidth = 1 / float(c.width);
	double inverseheight = 1 / float(c.height);
	float aspectratio = c.width / c.height;
	float angle = tan(3.14 * .5 * c.fov / 180.0f);
	glm::vec4 roish = c.toWorld * glm::vec4{ c.position, 0. };
	glm::vec3 ro = glm::vec3{ roish.x,roish.y,roish.z };
	int depth = 0;
	c.renderTarget.assign(c.width, c.height, 1, 3);
	double xx, yy;

	glm::mat4 transform = glm::translate(glm::mat4(1.0f), { .75,-.25,-4. });
	transform = glm::rotate(transform, glm::radians(-90.0f), { 1,0,0 });
	transform = glm::rotate(transform, glm::radians(-135.0f), { 0,0,1 });
	transform = glm::scale(transform, { .05,.05,.05 });

	mat1->matcolor = glm::vec3{ 0.77, 0.78, 0.78 };
	mat1->ref = 0;
	mat1->ior = 0;
	mat1->roughness = .2;
	mat1->metalness = 1.;

	std::string warning, error;
	bool loaded = tinyobj::LoadObj(mod1.attrib, mod1.shapes, mod1.materials, &warning, &error, "teapot.obj");

	mod1.transform = transform;
	mod1.invtransform = glm::inverse(transform);
	mod1.invtrastransform = glm::transpose(mod1.invtransform);
	mod1.modelnum = 0;
	mod1.modeltype = 0;
	mod1.mat = mat1;

	transform = glm::translate(glm::mat4(1.0f), { -.75,-.25,-4. });
	transform = glm::rotate(transform, glm::radians(-90.0f), { 1,0,0 });
	transform = glm::rotate(transform, glm::radians(45.0f), { 0,0,1 });
	transform = glm::scale(transform, { .05,.05,.05 });

	mat2->matcolor = glm::vec3{ 0.77, 0.78, 0.78 };
	mat2->ref = 0;
	mat2->roughness = .2;
	mat2->metalness = 0.;

	loaded = tinyobj::LoadObj(mod2.attrib, mod2.shapes, mod2.materials, &warning, &error, "teapot.obj");

	mod2.transform = transform;
	mod2.invtransform = glm::inverse(transform);
	mod2.invtrastransform = glm::transpose(mod2.invtransform);
	mod2.modelnum = 1;
	mod2.modeltype = 0;
	mod2.mat = mat2;

	mat3->matcolor = glm::vec3{ 0.,0.,1. };
	mat3->ref = 0.;
	mat3->ior = 1.5;
	mat3->gangle = 3.14 / 90;

	mod3.c0 = glm::vec3{ .75,-.25,-4. };
	mod3.sr = .75;
	mod3.modelnum = 2;
	mod3.modeltype = 1;
	mod3.mat = mat3;

	mat4->matcolor = glm::vec3{ 1.,1.,0. };
	mat4->ref = 0;
	mat4->gangle = 3.14 / 90;
	mat4->roughness = .5;
	mat4->metalness = 0;

	mod4.c0 = glm::vec3{ -100,-2,0 };
	mod4.c1 = glm::vec3{ 100,-1.5,-500 };
	mod4.modeltype = 2;
	mod4.modelnum = 0;
	mod4.mat = mat4;

	mat5->matcolor = glm::vec3{ 0.77, 0.78, 0.78 };
	mat5->ref = .3;
	mat5->ior = 0.;
	mat5->roughness = .2;
	mat5->metalness = 1.;

	mod5.c0 = glm::vec3{ -1.15,0,-6 };
	mod5.sr = 1;
	mod5.modelnum = 1;
	mod5.modeltype = 1;
	mod5.mat = mat5;

	mat6->matcolor = glm::vec3{ 0.77, 0.78, 0.78 };
	mat6->ref = 0.;
	mat6->ior = 0.;
	mat6->roughness = .2;
	mat6->metalness = 0.;

	mod6.c0 = glm::vec3{ 1.15, 0.,-6 };
	mod6.sr = 1.;
	mod6.modelnum = 2;
	mod6.modeltype = 1;
	mod6.mat = mat6;
	

	struct Node n1;
	struct Node n2;
	n1 = BVH(mod1, n1);
	n2 = BVH(mod2, n2);

	std::vector<struct Node> n;
	n.push_back(n1);
	n.push_back(n2);
	std::vector<struct model> mod;
	//mod.push_back(mod1);
	//mod.push_back(mod2);
	//mod.push_back(mod3);
	mod.push_back(mod4);
	mod.push_back(mod5);
	mod.push_back(mod6);

	std::srand(time(NULL));

	time_t start = time(0);
	for (int y = 0; y < c.height; y++)
	{
		for (int x = 0; x < c.width; x++)
		{
			//if (y == 150 && x == 150)
			//{
				//std::cout << "kill me" << std::endl;
				//Sleep(5000);
			//}
			
			glm::vec3 color = { 0.,0.,0. };
			for (int r = 0; r < s.rPerPixel; r++)
			{
				float xoffset = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
				float yoffset = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
				xoffset = xoffset - .5f;
				yoffset = yoffset - .5f;

				xx = (2 * ((x + 0.5f + xoffset) * inversewidth) - 1) * angle * aspectratio;
				yy = (1 - 2 * ((y + 0.5f + yoffset) * inverseheight)) * angle;
				glm::vec3 ray = { xx,yy,-1 };

				//Lens Math
				if (c.lens_radius > 0)
				{
					float fdistance = glm::distance(c.lookAt, c.position);
					float radius = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
					float angle = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)) * (2 * 3.1415926535);
					glm::vec2 plancoord = glm::vec2{ radius * cos(angle),radius * sin(angle) };
					glm::vec2 plens = c.lens_radius * plancoord;

					float ft = abs(fdistance / ray.z);
					glm::vec3 focusp = ray * ft;

					ro = glm::vec3{ plens.x, plens.y, 0 };
					ray = glm::normalize(focusp - ro);
				}

				glm::vec4 temp = c.toWorld * glm::vec4(ray, 0);
				ray = glm::normalize(glm::vec3(temp.x, temp.y, temp.z));

				
				color += Trace(ray, ro, depth, c, s, { 0., 0., 0. }, light, mod, n);
				
				//color = ray;

			}


			color = color / float(s.rPerPixel);
			buffer[x][y] = color;
		}
	}

	float max_luminance = 0;
	if (!niaveRien)
	{
		for (int y = 0; y < c.height; y++)
		{
			for (int x = 0; x < c.width; x++)
			{
				glm::vec3 color = buffer[x][y];
				float temp = glm::dot(color, glm::vec3(0.2126f, 0.7152f, 0.0722f));
				if (temp > max_luminance)
				{
					max_luminance = temp;
				}
			}
		}
	}

	//Tone Mapping
	for (int y = 0; y < c.height; y++)
	{
		for (int x = 0; x < c.width; x++)
		{
			glm::vec3 color = buffer[x][y];
			if (niaveRien)
			{
				color = color / (1.f + color);
			}
			else
			{
				float l_old = glm::dot(buffer[x][y], glm::vec3(0.2126f, 0.7152f, 0.0722f));
				float numerator = l_old * (1.0f + (l_old / (max_luminance * max_luminance)));
				float l_new = numerator / (1.0f + l_old);
				color = buffer[x][y] * (l_new / l_old);
			}

			color = glm::clamp(color, { 0.,0.,0. }, { 1.,1.,1. });

			c.renderTarget(x, y, 0, 0) = (unsigned char)(abs(color.r) * 255);
			c.renderTarget(x, y, 0, 1) = (unsigned char)(abs(color.g) * 255);
			c.renderTarget(x, y, 0, 2) = (unsigned char)(abs(color.b) * 255);
		}
	}

	time_t end = time(0);
	time_t delta = end - start;

	float hours = delta / (60.0f * 60.0f);
	int intHours = int(hours);

	float minutes = (hours - intHours) * 60.0f;
	int intMinutes = int(minutes);

	float seconds = (minutes - intMinutes) * 60.0f;
	int intSeconds = int(seconds);

	std::cout << "Render Time: " << intHours << ", " << intMinutes << ", " << intSeconds << std::endl;
	c.renderTarget.save("test.ppm");
}

