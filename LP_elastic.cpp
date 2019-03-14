//a simple linear spring network model for elastic deformation of soft material
// developed by Yu Zheng, modified by Yang Jiao
// updated date: 09/01/2018

#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>

const int DIM = 402;
const double L = 1;
const double SPACE = L / (DIM - 1);
const double PI = 3.1415926;
const double DIAGONAL = sqrt(SPACE * SPACE + SPACE * SPACE);
const double KS = 0.01; //spring coefficient
const double MOVELENGTH = 1;
const int STEPS = 5000;
// for circle shape
const double INITR = L * 0.2;
const double INCREASERATIO = 1 - SPACE / 2;
// for triangle shape, semicircle shape, star shape strange triangle shape
const double SIDE = L * 0.3;
const double INITRS = L * 0.2;
const double INITINCREASE = SPACE / 2; // initial increase ratio
const double FININCREASE = 0.3;		   // final increase ratio
const double DIS = 0.8;				   // the distance between center of circle and center of triangle

struct node
{
	double position[2];
	bool inside;
	bool inside2;
	double force[2];
};

void initialization(node nodes[DIM][DIM])
{
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			nodes[i][j].position[0] = (i - (DIM - 1) / 2.) * SPACE;
			nodes[i][j].position[1] = (j - (DIM - 1) / 2.) * SPACE;
		}
	}
}

double findDistance(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

//*********for circle shape*******************
void findCircle(node nodes[DIM][DIM])
{
	//find insider
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			double distance;
			distance = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], 0.0, 0.0);

			if (distance <= INITR)
			{
				nodes[i][j].inside = true;
			}
			else
				nodes[i][j].inside = false;
		}
	}
}
void growCircle(node nodes[DIM][DIM])
{
	//grow circle
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (nodes[i][j].inside)
			{
				nodes[i][j].position[0] = nodes[i][j].position[0] * INCREASERATIO;
				nodes[i][j].position[1] = nodes[i][j].position[1] * INCREASERATIO;
			}
		}
	}
}
//*********for circle shape*******************

//*********for triangle shape*****************
void findTriangle(node nodes[DIM][DIM])
{
	//find insider
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if ((nodes[i][j].position[1] >= y) && (nodes[i][j].position[1] <= (-slope * std::abs(nodes[i][j].position[0]) + interceptY)))
			{
				nodes[i][j].inside = true;
			}
			else
				nodes[i][j].inside = false;
		}
	}
}
void growTriangle(node nodes[DIM][DIM])
{
	//grow triangle
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	const double vertex1[2] = {0, interceptY};
	const double vertex2[2] = {-SIDE * 0.5, y};
	const double vertex3[2] = {SIDE * 0.5, y};

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (nodes[i][j].inside)
			{
				double distance1 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex1[0], vertex1[1]);
				double distance2 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex2[0], vertex2[1]);
				double distance3 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex3[0], vertex3[1]);
				double minDistance = distance1 > distance2 ? (distance2 > distance3 ? distance3 : distance2) : (distance1 > distance3 ? distance3 : distance1);

				nodes[i][j].position[0] = nodes[i][j].position[0] + nodes[i][j].position[0] * INITINCREASE * (minDistance / interceptY);
				nodes[i][j].position[1] = nodes[i][j].position[1] + nodes[i][j].position[1] * INITINCREASE * (minDistance / interceptY);
			}
		}
	}
}
//*********for triangle shape*****************

//*********for semicircle shape*****************
void findSemicircle(node nodes[DIM][DIM])
{
	//find insider
	int numInsiders = 0;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			double distance;
			distance = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], 0.0, -INITRS / 2.0);

			if (distance <= INITRS && nodes[i][j].position[1] >= -INITRS / 2.0)
			{
				nodes[i][j].inside = true;
				numInsiders++;
			}
			else
				nodes[i][j].inside = false;
		}
	}
}
void growSemicircle(node nodes[DIM][DIM])
{
	//grow semicircle
	const double vertex1[2] = {-INITRS, 0};
	const double vertex2[2] = {INITRS, 0};
	const double interceptY = INITRS;

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (nodes[i][j].inside)
			{
				double distance1 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex1[0], vertex1[1]);
				double distance2 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex2[0], vertex2[1]);
				double minDistance = distance1 > distance2 ? distance2 : distance1;

				nodes[i][j].position[0] = nodes[i][j].position[0] + nodes[i][j].position[0] * INITINCREASE * (minDistance / interceptY);
				nodes[i][j].position[1] = nodes[i][j].position[1] + nodes[i][j].position[1] * INITINCREASE * (minDistance / interceptY);
			}
		}
	}
}
//*********for semicircle shape*****************

//*********for star shape*****************
void findStar(node nodes[DIM][DIM])
{
	//find insider
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if ((nodes[i][j].position[1] >= y) && (nodes[i][j].position[1] <= (-slope * std::abs(nodes[i][j].position[0]) + interceptY)))
			{
				nodes[i][j].inside = true;
			}
			else
			{
				if ((nodes[i][j].position[1] <= -y) && (nodes[i][j].position[1] >= (slope * std::abs(nodes[i][j].position[0]) - interceptY)))
				{
					nodes[i][j].inside = true;
				}
				else
					nodes[i][j].inside = false;
			}
		}
	}
}
void growStar(node nodes[DIM][DIM])
{
	//grow triangle
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	const double vertexs[12][2] = {{0, interceptY},
								   {-SIDE * 0.5, y},
								   {SIDE * 0.5, y},
								   {0, -interceptY},
								   {-SIDE * 0.5, -y},
								   {SIDE * 0.5, -y},
								   {(-y - interceptY) / slope, -y},
								   {(y + interceptY) / slope, -y},
								   {(y + interceptY) / slope, y},
								   {(-y - interceptY) / slope, y},
								   {-interceptY / slope, 0},
								   {interceptY / slope, 0}};

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			double minDistance = L;
			if (nodes[i][j].inside)
			{
				for (int n = 0; n < 12; n++)
				{
					double distance = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertexs[n][0], vertexs[n][1]);
					if (distance < minDistance)
						minDistance = distance;
				}

				nodes[i][j].position[0] = nodes[i][j].position[0] + nodes[i][j].position[0] * INITINCREASE * minDistance;
				nodes[i][j].position[1] = nodes[i][j].position[1] + nodes[i][j].position[1] * INITINCREASE * minDistance;
			}
		}
	}
}
//*********for star shape*****************

//*********for strange triangle shape*****************
void findStriangle(node nodes[DIM][DIM])
{
	//find insider
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	const double circleCenter[3][2] = {{DIS * -0.8660, DIS * 0.5},
									   {DIS * 0.8660, DIS * 0.5},
									   {0, -DIS}};
	const double radius = findDistance(0, -DIS, SIDE * 0.5, y);

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			double distance1 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], circleCenter[0][0], circleCenter[0][1]);
			double distance2 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], circleCenter[1][0], circleCenter[1][1]);
			double distance3 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], circleCenter[2][0], circleCenter[2][1]);

			if ((nodes[i][j].position[1] >= y) && (nodes[i][j].position[1] <= (-slope * std::abs(nodes[i][j].position[0]) + interceptY)) && distance1 >= radius && distance2 >= radius && distance3 >= radius)
			{
				nodes[i][j].inside = true;
			}
			else
				nodes[i][j].inside = false;
		}
	}
}

void growStriangle(node nodes[DIM][DIM])
{
	//grow triangle
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double y = -SIDE * slope * (1.0 / 6.0);
	const double vertex1[2] = {0, interceptY};
	const double vertex2[2] = {-SIDE * 0.5, y};
	const double vertex3[2] = {SIDE * 0.5, y};

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (nodes[i][j].inside)
			{
				double distance1 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex1[0], vertex1[1]);
				double distance2 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex2[0], vertex2[1]);
				double distance3 = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], vertex3[0], vertex3[1]);
				double maxDistance = distance1 < distance2 ? (distance2 < distance3 ? distance3 : distance2) : (distance1 < distance3 ? distance3 : distance1);

				nodes[i][j].position[0] = nodes[i][j].position[0] - nodes[i][j].position[0] * INITINCREASE * (maxDistance / SIDE);
				nodes[i][j].position[1] = nodes[i][j].position[1] - nodes[i][j].position[1] * INITINCREASE * (maxDistance / SIDE);
			}
		}
	}
}
//*********for strange triangle shape*****************

void addcellTri(node nodes[DIM][DIM])
{
	const double slope = sqrt(3);
	const double interceptY = SIDE * slope * (2.0 / 6.0);
	const double interceptY2 = interceptY + 0.1;
	const double y1 = -SIDE * slope * (1.0 / 6.0);
	const double y2 = y1 - 0.05;
	const double vertex1[2] = {0, interceptY};
	const double vertex2[2] = {-SIDE * 0.5, y1};
	const double vertex3[2] = {SIDE * 0.5, y1};

	const double radius = SPACE * 4;
	const double density = 0.3;
	const double area = (interceptY2 - y2) * (interceptY2 - y2) / sqrt(3) - SIDE * (interceptY - y1) * 0.5;
	const int number = int(density * area / (PI * radius * radius));
	double positionCell1[number][2];

	int k1 = 0;
	srand((unsigned int)(time(NULL)));
	while (k1 < number)
	{
		double x = 0.4 * (std::rand() / double(RAND_MAX)) - 0.2;
		double y = 0.5 * (std::rand() / double(RAND_MAX)) - 0.2;
		if ((y <= (-slope * std::abs(x) + interceptY2)) && (y >= (-slope * std::abs(x) + interceptY)) || (y >= y2 && y <= y1))
		{
			bool isOverlap = false;
			int i = 0;

			while (!isOverlap && i < k1)
			{
				if (findDistance(x, y, positionCell1[i][0], positionCell1[i][1]) < 2 * radius)
					isOverlap = true;
				i++;
			}
			if (isOverlap == false)
			{
				positionCell1[k1][0] = x;
				positionCell1[k1][1] = y;
				k1++;
			}
		}
	}

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			for (int n = 0; n < number; n++)
			{
				double distance = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], positionCell1[n][0], positionCell1[n][1]);
				if ((distance >= radius - 0.5 * radius) && (distance <= radius + 0.5 * radius))
				{
					double changeX = 0.15 * (positionCell1[n][0] - nodes[i][j].position[0]);
					double changeY = 0.15 * (positionCell1[n][1] - nodes[i][j].position[1]);
					double rsquared = positionCell1[n][0] * positionCell1[n][0] + positionCell1[n][1] * positionCell1[n][1];
					double projectChangeX = positionCell1[n][0] * (changeX * positionCell1[n][0] + changeY * positionCell1[n][1]) / rsquared;
					double projectChangeY = positionCell1[n][1] * (changeX * positionCell1[n][0] + changeY * positionCell1[n][1]) / rsquared;

					nodes[i][j].position[0] = projectChangeX + nodes[i][j].position[0];
					nodes[i][j].position[1] = projectChangeY + nodes[i][j].position[1];
				}
			}
		}
	}

	std::ofstream fout2("test2.txt");
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!nodes[i][j].inside)
				fout2 << nodes[i][j].position[0] << " " << nodes[i][j].position[1] << std::endl;
		}
	}
	fout2.close();

	std::ofstream fout3("cellposition.txt");
	for (int i = 0; i < number; i++)
	{
		fout3 << positionCell1[i][0] << " " << positionCell1[i][1] << std::endl;
	}
	fout3.close();
}

void addcellCircle(node nodes[DIM][DIM])
{
	const double density = 0.1;
	const double radiusOutside = INITR + 0.1;
	const double radius = SPACE * 3;
	const double area = PI * radiusOutside * radiusOutside - PI * INITR * INITR;
	const int number = int(density * area / (PI * radius * radius));
	double positionCell1[number][2];

	int k1 = 0;
	srand((unsigned int)(time(NULL)));
	while (k1 < number)
	{
		double x = 0.6 * (std::rand() / double(RAND_MAX)) - 0.3;
		double y = 0.6 * (std::rand() / double(RAND_MAX)) - 0.3;
		if ((x * x + y * y <= radiusOutside * radiusOutside) && (x * x + y * y >= INITR * INITR))
		{
			bool isOverlap = false;
			int i = 0;

			while (!isOverlap && i < k1)
			{
				if (findDistance(x, y, positionCell1[i][0], positionCell1[i][1]) < 2 * radius)
					isOverlap = true;
				i++;
			}
			if (isOverlap == false)
			{
				positionCell1[k1][0] = x;
				positionCell1[k1][1] = y;
				k1++;
			}
		}
	}

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			for (int n = 0; n < number; n++)
			{
				double distance = findDistance(nodes[i][j].position[0], nodes[i][j].position[1], positionCell1[n][0], positionCell1[n][1]);
				if ((distance >= radius - 0.5 * radius) && (distance <= radius + 0.5 * radius))
				{
					double changeX = 0.01 * (positionCell1[n][0] - nodes[i][j].position[0]);
					double changeY = 0.01 * (positionCell1[n][1] - nodes[i][j].position[1]);
					double rsquared = positionCell1[n][0] * positionCell1[n][0] + positionCell1[n][1] * positionCell1[n][1];
					double projectChangeX = positionCell1[n][0] * (changeX * positionCell1[n][0] + changeY * positionCell1[n][1]) / rsquared;
					double projectChangeY = positionCell1[n][1] * (changeX * positionCell1[n][0] + changeY * positionCell1[n][1]) / rsquared;

					nodes[i][j].position[0] = projectChangeX + nodes[i][j].position[0];
					nodes[i][j].position[1] = projectChangeY + nodes[i][j].position[1];
				}
			}
		}
	}

	std::ofstream fout2("test2.txt");
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!nodes[i][j].inside)
				fout2 << nodes[i][j].position[0] << " " << nodes[i][j].position[1] << std::endl;
		}
	}
	fout2.close();

	std::ofstream fout3("cellposition.txt");
	for (int i = 0; i < number; i++)
	{
		fout3 << positionCell1[i][0] << " " << positionCell1[i][1] << std::endl;
	}
	fout3.close();
}

void findForce(node nodes[DIM][DIM])
{
	int i, j;
#pragma omp parallel for private(j)
	for (i = 1; i < DIM - 1; i++)
	{
		for (j = 1; j < DIM - 1; j++)
		{
			int indices[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

			nodes[i][j].force[0] = 0;
			nodes[i][j].force[1] = 0;
			if (!nodes[i][j].inside)
			{
				for (int k = 0; k < 4; k++)
				{
					double distance;
					double deltaX;
					double deltaY;
					double fxk, fyk;

					int di = indices[k][0], dj = indices[k][1];

					deltaX = nodes[i + di][j + dj].position[0] - nodes[i][j].position[0];
					deltaY = nodes[i + di][j + dj].position[1] - nodes[i][j].position[1];
					distance = sqrt(deltaX * deltaX + deltaY * deltaY);

					fxk = (deltaX / distance) * (distance - SPACE) * KS;
					fyk = (deltaY / distance) * (distance - SPACE) * KS;
					nodes[i][j].force[0] += fxk;
					nodes[i][j].force[1] += fyk;
				}
			}

			int indices2[4][2] = {{-1, 1}, {1, 1}, {-1, -1}, {1, -1}};
			if (!nodes[i][j].inside)
			{
				for (int k = 0; k < 4; k++)
				{
					double distance;
					double deltaX;
					double deltaY;
					double fxk, fyk;

					int di = indices2[k][0], dj = indices2[k][1];

					deltaX = nodes[i + di][j + dj].position[0] - nodes[i][j].position[0];
					deltaY = nodes[i + di][j + dj].position[1] - nodes[i][j].position[1];
					distance = sqrt(deltaX * deltaX + deltaY * deltaY);

					fxk = (deltaX / distance) * (distance - DIAGONAL) * KS;
					fyk = (deltaY / distance) * (distance - DIAGONAL) * KS;
					nodes[i][j].force[0] += fxk;
					nodes[i][j].force[1] += fyk;
				}
			}
		}
	}
}

void move(node nodes[DIM][DIM])
{
	int i, j;
#pragma omp parallel for private(j)
	for (i = 1; i < DIM - 1; i++)
	{
		for (j = 1; j < DIM - 1; j++)
		{
			if (!nodes[i][j].inside)
			{
				if (std::abs(MOVELENGTH * nodes[i][j].force[0]) > (0.05 * SPACE))
					nodes[i][j].position[0] += 0.05 * SPACE * nodes[i][j].force[0];
				else
					nodes[i][j].position[0] += MOVELENGTH * nodes[i][j].force[0];
				if (std::abs(MOVELENGTH * nodes[i][j].force[1]) > (0.05 * SPACE))
					nodes[i][j].position[1] += 0.05 * SPACE * nodes[i][j].force[1];
				else
					nodes[i][j].position[1] += MOVELENGTH * nodes[i][j].force[1];
			}
		}
	}
}

int main()
{
	double startTime, endTime;
	int stepsincrease = floor(log(1 + FININCREASE) / log(1 + INITINCREASE));
	startTime = omp_get_wtime();

	node nodes[DIM][DIM];
	initialization(nodes);

	findTriangle(nodes);

	std::ofstream fout1("Tri_originpositiontest.txt");
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!nodes[i][j].inside)
				fout1 << nodes[i][j].position[0] << " " << nodes[i][j].position[1] << std::endl;
		}
	}
	fout1.close();
	/*
	for (int nincrease = 0; nincrease < stepsincrease; nincrease++)
	{

		growStriangle(nodes);

		for (int nsteps = 0; nsteps < STEPS; nsteps++)
		{
			findForce(nodes);
			move(nodes);
		}
	}
*/
	addcellTri(nodes);
	/*
	findCircle(nodes);
	std::ofstream fout1("400circle_pull_originpositiontest.txt");
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!nodes[i][j].inside)
				fout1 << nodes[i][j].position[0] << " " << nodes[i][j].position[1] << std::endl;
		}
	}
	fout1.close();
	growCircle(nodes);
	for (int nsteps = 0; nsteps < STEPS; nsteps++)
	{
		findForce(nodes);
		move(nodes);
	}

	std::ofstream fout("400circle_pull_	positionstest.txt");
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			if (!nodes[i][j].inside)
				fout << nodes[i][j].position[0] << " " << nodes[i][j].position[1] << " " << nodes[i][j].force[0] << " " << nodes[i][j].force[1] << " " << sqrt(nodes[i][j].force[0] * nodes[i][j].force[0] + nodes[i][j].force[1] * nodes[i][j].force[1]) << std::endl;
		}
	}
	fout.close();
	*/
	endTime = omp_get_wtime();
	std::cout << "Total time: " << endTime - startTime << " s" << std::endl;

	return 0;
}
