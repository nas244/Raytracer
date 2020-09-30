#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/transform.hpp"
#include "CImg.h"
#define TINYOBJLOADER_IMPLEMENTATION
//#include "tiny_obj_loader.h"
#include <vector>

namespace tinyobj
{
	struct attrib_t;
	struct shape_t;
	struct material_t;
}

struct camera
{
	int width = 1920;
	int height = 1080;
	float fov = 45;
	glm::vec3 position = { 0.,0.,0. };
	glm::vec3 lookAt = { 0.,0.,-5. };
	glm::vec3 up = { 0.,1.,0. };
	glm::mat4 toWorld = glm::lookAt(position, lookAt, up);
	cimg_library::CImg<unsigned char> renderTarget;

};

struct directionalLight
{
	glm::vec3 direction = { .5,-1.,-1. };
	glm::vec3 color = { 1.,1.,1. };
};

struct scene
{
	int max_depth = 2.;
	int rPerPixel = 1;
};

struct materials
{
	glm::vec3 matcolor;
	float ref = 0.;
	glm::vec3 emissive;
};

struct model
{
	std::vector<tinyobj::shape_t>* shapes = new std::vector<tinyobj::shape_t>();
	std::vector<tinyobj::material_t>* materials = new std::vector<tinyobj::material_t>();
	//std::vector<tinyobj::shape_t>* attrib = new tinyobj::attrib_t();
};

struct hitInfo
{
	glm::vec3 hitPosition;
	glm::vec3 hitNormal;
	struct materials* hitMaterial;
	float hitDistance;
};


// Ray triangle intersection using the Moller-Trumbore algorithm (scratchapixel.com)
bool TestTriIntersection(glm::vec3 ro, glm::vec3 ray, struct hitInfo hit)
{
	glm::vec3 v0 = { 5.,0.,-5. };
	glm::vec3 v1 = { -5.,0.,-5. };
	glm::vec3 v2 = { 0.,5.,-5. };

	glm::vec3 v0v1 = v1 - v0;
	glm::vec3 v0v2 = v2 - v0;
	glm::vec3 pvec = glm::cross(ray, v0v2);
	float det = glm::dot(v0v1, pvec);

	if (fabs(det) < 0.0001) return false;

	float invDet = 1.0f / det;

	glm::vec3 tvec = ro - v0;
	float u = glm::dot(tvec, pvec) * invDet;
	if (u < 0 || u > 1) return false;

	glm::vec3 qvec = glm::cross(tvec, v0v1);
	float v = glm::dot(ray, qvec) * invDet;

	if (v < 0 || u + v > 1) return false;

	hit.hitDistance = glm::dot(v0v2, qvec) * invDet;

	if (hit.hitDistance < 0.0001) return false;

	//hitInfo.hitPosition = rayOrigin + rayDirection * hitInfo.hitDistance;
	//hitInfo.hitMaterial = m_material;

	//h.hitNormal = glm::normalize((1.0f - u - v) * n0 + u * n1 + v * n2);

	hit.hitNormal = glm::vec3{ 0.,0.,-1. };
	//glm::vec4 tempNormal = m_transform * glm::vec4(interpolatedNormal, 0);
	//hitInfo.hitNormal = { tempNormal.x, tempNormal.y, tempNormal.z };

	return true;
}

glm::vec3 Trace(glm::vec3 ray, glm::vec3 ro, int depth, struct scene s, glm::vec3 color, struct directionalLight light)
{
	if (depth >= s.max_depth)
	{
		return color;
	}

	struct hitInfo hit;

	bool intersectsObject = TestTriIntersection(ro, ray, hit);

	if (!intersectsObject)
	{
		return color;
	}

	struct hitInfo occulsion;

	bool occluded = TestTriIntersection(hit.hitPosition + 0.001f * light.direction, light.direction, occulsion);

	if (occluded)
	{
		return color;
	}

	color = glm::vec3{ 1.,0.,0. };

	return color;
}

int main(int argc, char* argv[])
{
	struct camera c;
	struct directionalLight light;
	struct scene s;

	float inversewidth = 1 / c.width;
	float inverseheight = 1 / c.height;
	float aspectratio = c.width / c.height;
	float angle = tan(3.14 * .5 * c.fov / 180.0f);
	glm::vec4 ro = c.toWorld * glm::vec4{ c.position, 0.};
	int depth = 0;
	c.renderTarget.assign(c.width, c.height, 1, 3);

	glm::mat4 transform = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -10));
	transform = glm::rotate(transform, glm::radians(-90.0f), { 1,0,0 });

	//std::string warning, error;

	//bool loaded = tinyobj::LoadObj(attrib, shapes, materials, &warning, &error, file);

	for (int y = 0; y < c.height; y++)
	{
		for (int x = 0; x < c.width; x++)
		{
			glm::vec3 color = { 0.,0.,0. };
			for (int r = 0; r < s.rPerPixel; r++)
			{
				float xx = (2 * ((x + 0.5f) * inversewidth) - 1) * angle * aspectratio;
				float yy = (1 - 2 * ((y + 0.5f) * inverseheight)) * angle;
				glm::vec3 ray = { xx,yy,-1 };

				glm::vec4 temp = glm::normalize(c.toWorld * glm::vec4(ray, 1));
				ray = glm::normalize(glm::vec3(temp.x, temp.y, temp.z));

				color += Trace(ray, { ro.x, ro.y, ro.z }, depth, s, {0., 0., 0.}, light);
			}

			color = color / float(s.rPerPixel);

			color = glm::clamp(color, { 0.,0.,0. }, { 1.,1.,1. });

			c.renderTarget(x, y, 1, 0) = (unsigned char)(color.r * 255);
			c.renderTarget(x, y, 1, 1) = (unsigned char)(color.g * 255);
			c.renderTarget(x, y, 1, 2) = (unsigned char)(color.b * 255);

			c.renderTarget.save("test");
		}
	}
}