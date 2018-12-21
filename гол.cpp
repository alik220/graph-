#include <vector>
#include <cmath>
#include <stack>
#include <algorithm>
#include <iostream>
#include <map>
#include <cstdlib>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

/*
//**********************************************РИСУЕМ




const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);

using namespace std;

const int width  = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
{
	bool steep = false;

	if (abs(x0 - x1) < abs(y0 - y1))
	{
		swap(x0, y0);
		swap(x1, y1);
		steep = true;
	}
	if (x0 > x1)
	{
		swap(x0, x1);
		swap(y0, y1);
	}

	int dx = x1 - x0;
	int dy = y1 - y0;
	int derror2 = abs(dy) * 2;
	int error2 = 0;
	int y = y0;

	for (int x = x0; x <= x1; x++)
	{
		if (steep)
		{
			image.set(y, x, color);
		}
		else
		{
			image.set(x, y, color);
		}

		error2 += derror2;

		if (error2 > dx)
		{
			y += (y1 > y0 ? 1 : -1);
			error2 -= dx * 2;
		}
	}
}

void triangle(Vec2i p1, Vec2i p2, Vec2i p3, TGAImage& image, TGAColor color)
{
	line(p1.x, p1.y, p2.x, p2.y, image, color);
	line(p2.x, p2.y, p3.x, p3.y, image, color);
	line(p3.x, p3.y, p1.x, p1.y, image, color);
}

int main()
{
	Model* model = new Model("african_head.obj");

	TGAImage image(width, height, TGAImage::RGB);

// v - точки
// f - грани. if возьмем 1ое число из каждой строчки - треуг-ки
		for (int i = 0; i < model->nfaces(); i++)
		{
			std::vector<int> face = model->face(i);

			for (int j = 0; j < 3; j++)
			{
				Vec3f v0 = model->vert(face[j]);
				Vec3f v1 = model->vert(face[(j + 1) % 3]);

				int x0 = (v0.x + 1.)*width / 2.;
				int y0 = (v0.y + 1.)*height / 2.;
				int x1 = (v1.x + 1.)*width / 2.;
				int y1 = (v1.y + 1.)*height / 2.;
				line(x0, y0, x1, y1, image, white);

				//cout << x0 << " " << y0 << " ; " << x1 << " " << y1 << endl;

			}
		}

		image.flip_vertically();
		image.write_tga_file("fil.tga");

		delete model;
		return 0;
	}




*/


//********************************************************КРАСИМ
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);

const int width = 2000;
const int height = 2000;

TGAImage image(width, height, TGAImage::RGBA);
Model *model = NULL;

using namespace std;




//t0, t1, t2 - т., отсортированные по возрастанию координаты у

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)

{
	if (t0.y == t1.y && t0.y == t2.y)
		return;

	if (t0.y > t1.y)
		swap(t0, t1);

	if (t0.y > t2.y)
		swap(t0, t2);

	if (t1.y > t2.y)
		swap(t1, t2);

	int total_height = t2.y - t0.y;

	for (int i = 0; i < total_height; i++)
	{
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;

		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;

		float alpha = (float)i / total_height;

		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;

		Vec2i A = t0 + (t2 - t0)*alpha;

		Vec2i B = second_half ? t1 + (t2 - t1)*beta : t0 + (t1 - t0)*beta;

		if (A.x > B.x)
			swap(A, B);

		for (int j = A.x; j <= B.x; j++)
		{
			image.set(j, t0.y + i, color);
		}
	}
}






int main(int argc, char **argv) {



	Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
	Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
	Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };

	triangle(t0[0], t0[1], t0[2], image, red);
	triangle(t1[0], t1[1], t1[2], image, white);
	triangle(t2[0], t2[1], t2[2], image, green);

	   	   

	//pаскраска 

	if (2 == argc)
	{
		model = new Model(argv[1]);
	}

	else
	{
		model = new Model("african_head.obj");
	}

	TGAImage image(width, height, TGAImage::RGB);


	for (int i = 0; i < model->nfaces(); i++)
	{
		vector<int> face = model->face(i);

		Vec2i screen_coords[3];

		for (int j = 0; j < 3; j++)
		{
			Vec3f world_coords = model->vert(face[j]);
			screen_coords[j] = Vec2i((world_coords.x + 1.)*width / 2., (world_coords.y + 1.)*height / 2.);
		}

		triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
	}


	image.flip_vertically();
	image.write_tga_file("output.tga");


	return 0;



}














/*void Bresenham_algorithm(int x0, int y0, int x1, int y1, TGAColor color) {
	//TGAImage image(100, 100, TGAImage::RGB);
	bool linedown = false;
	bool goodrelay = true;
	int deltax = abs(x1 - x0);
	int deltay = abs(y1 - y0);
	int I = 2;
	double m = (I * deltay) / (double) deltax;
	double w = I - m;
	double e = I / 2;
	if (x0 > x1) {
		swap(x0, x1);
		swap(y0, y1);
	}
	int maxY = y0 + y1;
	if (y1 < y0 & x0 < x1) {
		swap(y1, y0);
		linedown = true;
	}
	if (deltay > deltax) {
		swap(deltay, deltax);
		swap(y1, x1);
		swap(y0, x0);
		goodrelay = false;
	}
	int error = 0;
	int deltaerr = deltay;
	int x = x0;
	int y = y0;
	int diry = y1 - y0;
	if (diry > 0)
		diry = 1;
	if (diry < 0)
		diry = -1;
	//for (int x = x0; x <= x1; x++) {
	while (x < x1) {
		if (!goodrelay)
			swap(y, x);
		if (linedown) {
			image.set(x, maxY - y, TGAColor(255, 255, 255, e * 255));
		}
		else {
			image.set(x, y, TGAColor(255, 255, 255, e * 255));
		}
		if (!goodrelay)
			swap(y, x);
		error = error + deltaerr;
		if (2 * error >= deltax) {
			y = y + diry;
			error = error - deltax;
			e -= w;
			x++;
		}
		else {
			e += m;
			x++;
		}
		//cout << "e = " << e << endl ;
	}
	//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//image.write_tga_file("output.tga");
}
*/








/*
 //Проволочны рендер
	if (2==argc) {
		model = new Model(argv[1]);
	} else {
		model = new Model("african_head.obj");
	}

	for (int i=0; i<model->nfaces(); i++) {
		std::vector<int> face = model->face(i);
		for (int j=0; j<3; j++) {
			Vec3f v0 = model->vert(face[j]);
			Vec3f v1 = model->vert(face[(j+1)%3]);
			int x0 = (v0.x+1.)*width/2.;
			int y0 = (v0.y+1.)*height/2.;
			int x1 = (v1.x+1.)*width/2.;
			int y1 = (v1.y+1.)*height/2.;
			//cout << x0 << " " << y0 << " ; " << x1 << " " << y1 << endl;
			drawLineBresenham(x0, y0, x1, y1, white);
		}
	}


*/




/*vector<vector<double>> finalCoords;


for (int i = 0; i < model->nfaces(); i++)
{
	vector<int> face = model->face(i);

	for (int j = 0; j < 3; j++)
	{
		Vec3f coords = model->vert(face[j]);
		finalCoords.push_back({coords.x, coords.y, coords.z});
	}
}

for (int i = 0; i < finalCoords.size(); i += 3)
{
	vector<double> p1 = finalCoords[i];
	vector<double> p2 = finalCoords[i + 1];
	vector<double> p3 = finalCoords[i + 2];

	Vec2i pp1 = Vec2i((p1[0] + 1.0) * width / 2.0, (p1[1] + 1.0) * height / 2.0);
	Vec2i pp2 = Vec2i((p2[0] + 1.0) * width / 2.0, (p2[1] + 1.0) * height / 2.0);
	Vec2i pp3 = Vec2i((p3[0] + 1.0) * width / 2.0, (p3[1] + 1.0) * height / 2.0);

	triangle(pp1, pp2, pp3, image, white);
}
*/
