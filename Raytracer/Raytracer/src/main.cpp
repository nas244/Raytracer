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
	int width = 300;
	int height = 300;
	float fov = 45;
	glm::vec3 position = { 0.,0.,50. };
	glm::vec3 lookAt = { 0.,0.,-5. };
	glm::vec3 up = { 0.,1.,0. };
	glm::mat4 toWorld = glm::lookAt(position, lookAt, up);
	cimg_library::CImg<unsigned char> renderTarget;

};

struct directionalLight
{
	glm::vec3 direction = { 0.,1.,1. };
	glm::vec3 color = { 1.,1.,1. };
	float intensity = 1;
};

struct scene
{
	int max_depth = 2.;
	int rPerPixel = 1;
};

struct material
{
	glm::vec3 matcolor;
	float ref = 0.;
	glm::vec3 emissive;
};

struct model
{
	std::vector<tinyobj::shape_t>* shapes = new std::vector<tinyobj::shape_t>();
	std::vector<tinyobj::material_t>* materials = new std::vector<tinyobj::material_t>();
	tinyobj::attrib_t* attrib = new tinyobj::attrib_t();
};

struct hitInfo
{
	glm::vec3 hitPosition;
	glm::vec3 hitNormal;
	struct materials* hitMaterial;
	float hitDistance = 1000000;
	float hitMisses;
	bool didHit;
	float u;
	float v;
	float w;
	int index;
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
struct hitInfo TestBoxIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo hit, glm::vec3 c1, glm::vec3 c2)
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
		hit.didHit = false;
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
		hit.didHit = false;
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

	hit.hitDistance = t_min;
	hit.didHit = true;
	hit.hitPosition = glm::vec3{0.,0.,0.} + ray * hit.hitDistance;
	return hit;
}

struct hitInfo TestSphIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo hit)
{
	float tca = 0.;

	glm::vec3 so = { 0.,0.,-5 };
	int sr = 1;
	//sphere intersect
	glm::vec3 l = so - ro;
	tca = glm::dot(l, ray);
	if (tca < 0.)
	{
		hit.didHit = false;
		return hit;
	}
	float d2 = glm::dot(l, l) - tca * tca;
	if (d2 > sr * sr)
	{
		hit.didHit = false;
		return hit;
	}
	float thc = sqrt(sr * sr - d2);
	if (tca - thc < tca + thc && tca - thc > 0)
	{
		hit.hitDistance = tca - thc;
	}
	else if (tca - thc > tca + thc && tca + thc > 0)
	{
		hit.hitDistance = tca + thc;
	}

	hit.hitPosition = ro + ray * hit.hitDistance;

	hit.hitNormal = so - hit.hitPosition;

	hit.hitNormal = glm::normalize(hit.hitNormal);
	hit.didHit = true;

	return hit;
}


// Ray triangle intersection using the Moller-Trumbore algorithm (scratchapixel.com)
struct hitInfo TestTriIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo hit, glm::vec3 v0, glm::vec3 v1, glm::vec3 v2 )
{
	glm::vec3 v0v1 = v1 - v0;
	glm::vec3 v0v2 = v2 - v0;
	glm::vec3 pvec = glm::cross(ray, v0v2);
	float det = glm::dot(v0v1, pvec);

	if (fabs(det) < 0.0001)
	{
		hit.didHit = false;
		return hit;
	}

	float invDet = 1.0f / det;

	glm::vec3 tvec = ro - v0;
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0 || u > 1)
	{
		hit.didHit = false;
		return hit;
	}

	glm::vec3 qvec = glm::cross(tvec, v0v1);
	float v = glm::dot(ray, qvec) * invDet;

	if (v < 0 || u + v > 1)
	{
		hit.didHit = false;
		return hit;
	}

	hit.hitDistance = glm::dot(v0v2, qvec) * invDet;

	if (hit.hitDistance < 0.0001)
	{
		hit.didHit = false;
		return hit;
	}
		
	//hitInfo.hitPosition = rayOrigin + rayDirection * hitInfo.hitDistance;
	//hitInfo.hitMaterial = m_material;

	//h.hitNormal = glm::normalize((1.0f - u - v) * n0 + u * n1 + v * n2);

	hit.didHit = true;
	hit.u = u;
	hit.v = v;
	hit.hitPosition = glm::vec3{ 0.,0.,0. } + ray * hit.hitDistance;

	//hit.hitNormal = glm::vec3{ 0.,0.,-1. };
	//glm::vec4 tempNormal = m_transform * glm::vec4(interpolatedNormal, 0);
	//hitInfo.hitNormal = { tempNormal.x, tempNormal.y, tempNormal.z };

	return hit;
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct hitInfo traverse(struct model mod, struct Node n, struct hitInfo hit, glm::vec3 ro, glm::vec3 ray)
{
	struct hitInfo box;
	box = TestBoxIntersection(ro, ray, hit, n.c0, n.c1);

	if (!box.didHit)
	{
		hit.didHit = box.didHit;
		return hit;
	}

	struct hitInfo triangle;
	if (n.level == n.maxLevel)
	{
		glm::vec3 v0, v1, v2;
		int closest = -1;
		float closestHit = 1000000;
		for (int j = 0; j < n.index->size(); j += 1)
		{
			int index_0 = (*n.indices)[(*n.index)[j]].vertex_index * 3;
			v0.x = mod.attrib->vertices[index_0];
			v0.y = mod.attrib->vertices[index_0 + 1];
			v0.z = mod.attrib->vertices[index_0 + 2];
			int normalIndex_0 = (*n.indices)[(*n.index)[j]].normal_index * 3;
			glm::vec3 normal_0 = { mod.attrib->normals[normalIndex_0],mod.attrib->normals[normalIndex_0 + 1],mod.attrib->normals[normalIndex_0 + 2] };

			int index_1 = (*n.indices)[(*n.index)[j] + 1].vertex_index * 3;
			v1.x = mod.attrib->vertices[index_1];
			v1.y = mod.attrib->vertices[index_1 + 1];
			v1.z = mod.attrib->vertices[index_1 + 2];
			int normalIndex_1 = (*n.indices)[(*n.index)[j] + 1].normal_index * 3;
			glm::vec3 normal_1 = { mod.attrib->normals[normalIndex_1],mod.attrib->normals[normalIndex_1 + 1],mod.attrib->normals[normalIndex_1 + 2] };

			int index_2 = (*n.indices)[(*n.index)[j] + 2].vertex_index * 3;
			v2.x = mod.attrib->vertices[index_2];
			v2.y = mod.attrib->vertices[index_2 + 1];
			v2.z = mod.attrib->vertices[index_2 + 2];
			int normalIndex_2 = (*n.indices)[(*n.index)[j] + 2].normal_index * 3;
			glm::vec3 normal_2 = { mod.attrib->normals[normalIndex_2],mod.attrib->normals[normalIndex_2 + 1],mod.attrib->normals[normalIndex_2 + 2] };

			triangle = TestTriIntersection(ro, ray, hit, v0, v1, v2);

			if (triangle.didHit)
			{
				if (triangle.hitDistance < hit.hitDistance)
				{
					hit.index = j;
					hit.hitPosition = triangle.hitPosition;
					hit.hitDistance = triangle.hitDistance;
					hit.u = triangle.u;
					hit.v = triangle.v;
					hit.w = 1 - (hit.u + hit.v);
					hit.hitNormal = hit.w * normal_0 + hit.u * normal_1 + hit.v * normal_2;
					hit.hitNormal = glm::normalize(hit.hitNormal);
				}
			}
		}		

		return hit;
	}

	hit = traverse(mod, (*n.child1), hit, ro, ray);
	hit = traverse(mod, (*n.child2), hit, ro, ray);
	//hit = traverse(mod, (*n.child3), hit, ro, ray);
	//hit = traverse(mod, (*n.child4), hit, ro, ray);
	//hit = traverse(mod, (*n.child5), hit, ro, ray);
	//hit = traverse(mod, (*n.child6), hit, ro, ray); 
	//hit = traverse(mod, (*n.child7), hit, ro, ray);
	//hit = traverse(mod, (*n.child8), hit, ro, ray);

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

		if ((v0.x <= max_x && v0.x >= min_x) && (v0.y <= max_y && v0.y >= min_y) && (v0.z <= max_z && v0.z >= min_z))
		{
			withinbounds = true;
		}
		if ((v1.x <= max_x && v1.x >= min_x) && (v1.y <= max_y && v1.y >= min_y) && (v1.z <= max_z && v1.z >= min_z))
		{
			withinbounds = true;
		}
		if ((v2.x <= max_x && v2.x >= min_x) && (v2.y <= max_y && v2.y >= min_y) && (v2.z <= max_z && v2.z >= min_z))
		{
			withinbounds = true;
		}
		
		if (withinbounds)
		{
			index->push_back(j);
		}
	}
	
	n1->index = index;

	float average_x = (max_x + min_x) / 2;
	float average_y = (max_y + min_y) / 2;
	float average_z = (max_z + min_z) / 2;

	n1->child1 = subdivide(mod, (*n1), glm::vec3{ max_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	n1->child2 = subdivide(mod, (*n1), glm::vec3{ max_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child3 = subdivide(mod, (*n1), glm::vec3{ max_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child4 = subdivide(mod, (*n1), glm::vec3{ max_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child5 = subdivide(mod, (*n1), glm::vec3{ min_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child6 = subdivide(mod, (*n1), glm::vec3{ min_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child7 = subdivide(mod, (*n1), glm::vec3{ min_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n1->child8 = subdivide(mod, (*n1), glm::vec3{ min_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });


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

	n.child1 = subdivide(mod, n, glm::vec3{ max_x+1, max_y+1, max_z+1 }, glm::vec3{ min_x-1, average_y, min_z-1 });
	n.child2 = subdivide(mod, n, glm::vec3{ min_x-1, min_y-1, min_z-1 }, glm::vec3{ max_x+1, average_y, max_z+1 });
	//n.child3 = subdivide(mod, n, glm::vec3{ max_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n.child4 = subdivide(mod, n, glm::vec3{ max_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	//n.child5 = subdivide(mod, n, glm::vec3{ min_x, max_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n.child6 = subdivide(mod, n, glm::vec3{ min_x, max_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	//n.child7 = subdivide(mod, n, glm::vec3{ min_x, min_y, max_z }, glm::vec3{ average_x, average_y, average_z });
	//n.child8 = subdivide(mod, n, glm::vec3{ min_x, min_y, min_z }, glm::vec3{ average_x, average_y, average_z });
	

	return n;
}

glm::vec3 Trace(glm::vec3 ray, glm::vec3 ro, int depth, struct scene s, glm::vec3 color, struct directionalLight light, struct hitInfo hit, struct model mod, struct Node n)
{
	if (depth >= s.max_depth)
	{
		return color;
	}


	hit.hitMisses = 0;


	glm::vec3 v0 = { 0.,0.,0. };
	glm::vec3 v1 = { 0.,0.,0. };
	glm::vec3 v2 = { 0.,0.,0. };
	int closest = -1;
	float closestHit = 1000000;

	hit = traverse(mod, n, hit, ro, ray);

	
	if (hit.hitDistance == 1000000)
	{
		return glm::vec3{ 1.,1.,1. };
	}

	



	
	

	//hit = TestSphIntersection(ro, ray, hit);

	//if (!hit.didHit)
	//{
		//hit.hitMisses += 1;
		//return glm::vec3{ 1.,1.,1. };
	//}

	//struct hitInfo occulsion;

	//occulsion = TestSphIntersection(hit.hitPosition + 0.001f * light.direction, light.direction, occulsion);

	//if (!occulsion.didHit)
	//{
		//return color;
	//}

	light.direction = glm::normalize(light.direction);

	float diff = glm::dot(hit.hitNormal, light.direction);

	if (diff < 0) 
	{
		diff = 0;
	}

	glm::vec3 diffuse = diff * light.intensity * light.color;
	color = glm::vec3{ 1.,0.,0. } * diffuse;

	//color = glm::vec3{ 1.,0.,0. };

	return color;
}




int main(int argc, char* argv[])
{
	struct camera c;
	struct directionalLight light;
	struct scene s;
	struct hitInfo hit;
	struct model mod;

	double inversewidth = 1 / float(c.width);
	double inverseheight = 1 / float(c.height);
	float aspectratio = c.width / c.height;
	float angle = tan(3.14 * .5 * c.fov / 180.0f);
	glm::vec4 ro = c.toWorld * glm::vec4{ c.position, 0. };
	int depth = 0;
	c.renderTarget.assign(c.width, c.height, 1, 3);
	double xx, yy;

	glm::mat4 transform = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -10));
	transform = glm::rotate(transform, glm::radians(-90.0f), { 1,0,0 });

	std::string warning, error;
	bool loaded = tinyobj::LoadObj(mod.attrib, mod.shapes, mod.materials, &warning, &error, "teapot.obj");

	struct Node n;
	n = BVH(mod, n);


	for (int y = 0; y < c.height; y++)
	{
		for (int x = 0; x < c.width; x++)
		{
			//if (y == 150 && x == 149)
			//{
				//std::cout << "kill me" << std::endl;
				//Sleep(5000);
			//}
			glm::vec3 color = { 0.,0.,0. };
			for (int r = 0; r < s.rPerPixel; r++)
			{

				xx = (2 * ((x + 0.5f) * inversewidth) - 1) * angle * aspectratio;
				yy = (1 - 2 * ((y + 0.5f) * inverseheight)) * angle;
				glm::vec3 ray = { xx,yy,-1 };

				//glm::vec4 temp = glm::normalize(c.toWorld * glm::vec4(ray, 1));
				//ray = glm::normalize(glm::vec3(temp.x, temp.y, temp.z));
				ray = glm::normalize(ray);

				color += Trace(ray, { ro.x, ro.y, ro.z }, depth, s, { 0., 0., 0. }, light, hit, mod, n);
				//color = ray;

			}

			color = color / float(s.rPerPixel);

			color = glm::clamp(color, { 0.,0.,0. }, { 1.,1.,1. });



			c.renderTarget(x, y, 0, 0) = (unsigned char)(abs(color.r) * 255);
			c.renderTarget(x, y, 0, 1) = (unsigned char)(abs(color.g) * 255);
			c.renderTarget(x, y, 0, 2) = (unsigned char)(abs(color.b) * 255);

			//std::cout << "x: " << x << " y, " << y << " misses: " << color.r << std::endl;




		}
	}
	c.renderTarget.save("test.ppm");
}

