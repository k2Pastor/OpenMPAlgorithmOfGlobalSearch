
#define _USE_MATH_DEFINES
#include <iostream>
#include <queue>
#include <list>
#include <iterator>
#include <cmath>
#include <ctime>
#include <omp.h>

using namespace std;

struct TestPoint
{
	double x;
	double z;
	TestPoint(double _x = 0, double _z = 0) : x(_x), z(_z) { }
};

struct Interval
{
	double R1;
	TestPoint* pLeft;
	TestPoint* pRight;
	Interval(double _R = 0, TestPoint *tpLeft = NULL, TestPoint *tpRight = NULL)
		: R1(_R), pLeft(tpLeft), pRight(tpRight) { }

	Interval (const Interval& i) {
		R1 = i.R1;
		pLeft = i.pLeft;
		pRight = i.pRight;
	}

	Interval& operator=(const Interval& i) {
		R1 = i.R1;
		pLeft = i.pLeft;
		pRight = i.pRight;
		return *this;

	}


};

bool operator<(const Interval& i1, const Interval& i2) { return (i1.R1 < i2.R1) ? true : false; }

double f(double x)
{
	//return 2 * (x - 3.0) * (x - 3.0) + exp(x * x / 2.0); // совпадение;
	return sin(x) + sin((10 * x) / 3); // совпадение;
	//return x * x;
	//return ((3.0 * x - 1.4) * sin(18.0 * x)); // совпадение;
	//return -1.0 * (x + sin(x)) * exp(-1.0 * x * x); // совпадение;
	//return sin(x) + sin((10 * x) / 3) + log(x) - 0.84 * x + 3; // совпадение;
	//return -1.0 * sin(2 * M_PI * x) * exp(-1.0 * x); // почти совпадение;
	//return (x * x - 5 * x + 6) / (x * x + 1); // совпадение;
	//return -x + sin(3 * x) - 1; // совпадение;
}

/*double f(double x)
{
double sum = 0.0;
for (int k = 1; k <= 5; k++)
{
sum += k * sin((k + 1) * x + k);
}

return -1.0 * sum;
} */

double ComputeR(const TestPoint& tpLeft, const TestPoint& tpRight, double _m)
{
	double diffX = tpRight.x - tpLeft.x;
	double diffZ = tpRight.z - tpLeft.z;
	return _m * diffX + diffZ * diffZ / (_m * diffX) - 2 * (tpRight.z + tpLeft.z);
}

TestPoint* InsertUp_List(list<TestPoint> *ltp, TestPoint *tpk)
{
	list<TestPoint>::iterator itLeft, itRight;
	itLeft = itRight = (*ltp).begin();

	while ((itRight != (*ltp).end()) && (itRight->x < (*tpk).x))
	{
		itLeft = itRight;
		itRight++;
	}

	(*ltp).insert(itRight, (*tpk));
	itLeft++;

	return &(*itLeft);
}


int main()
{
	list<TestPoint> testPoints; // точки испытаний;
	list<TestPoint>::iterator itLeft, itRight;
	priority_queue<Interval> Queue;
	double a, b; // границы интервала;
	TestPoint* tpt;
	Interval CharacteristicInterval[4];
	TestPoint Points[4];
	double accuracy; // точность алгоритма;
	int iterations; // количество итераций;
	int k = 0; // количество испытаний;
	double M;
	double m = -1.0; // константа Липшица;
	double r = 2.0; // заданный параметр метода;
	double GlobalMin, DotOfGlobalMin;
	double timeStart, timeEnd;


	cout << "Enter the ends of the segment [a;b] : " << endl;
	cin >> a >> b;
	cout << "Enter number of the iterations: " << endl;
	cin >> iterations;
	cout << "Enter value of the accuracy: " << endl;
	cin >> accuracy;

	timeStart = clock();
	double pr = (b - a) / 4;
	for (int i = 0; i < 4; i++) {
		testPoints.push_back(TestPoint(a + pr * i, f(a + pr * i)));
	}

	testPoints.push_back(TestPoint(b, f(b)));



	//tp1.z = f(tp1.x);
	//tp2.z = f(tp2.x);

	//testPoints.push_back(tp1);
	//testPoints.push_back(tp2);
	k = 0;



	do
	{
		M = 0.0;
		itLeft = itRight = testPoints.begin();
		++itRight;

		double Old_m = m;

		while (itRight != testPoints.end())
		{
			double max_M = fabs((itRight->z - itLeft->z) / (itRight->x - itLeft->x));
			if (M < max_M)
			{
				M = max_M;
			}

			itRight++;
			itLeft++;
		}

		if (M > 0.0)
		{
			m = r * M;
		}
		else {
			m = 1.0;
		}

		if (Old_m != m)
		{
			Queue = priority_queue<Interval>();

			itLeft = itRight = testPoints.begin();
			itRight++;

			while (itRight != testPoints.end())
			{
				Queue.push(Interval(ComputeR(*itLeft, *itRight, m), &(*itLeft), &(*itRight)));
				++itLeft;
				++itRight;
			}

		}

		/*while (itRight != testPoints.end())
		{
			Queue.push(Interval(ComputeR(*itLeft, *itRight, m), &(*itLeft), &(*itRight)));
			++itLeft;
			++itRight;
		}

		if (itRight == testPoints.end())
		{
			(itRight--);
		}

		CharacteristicInterval = Queue.top();
		Queue.pop(); */

		CharacteristicInterval[0] = Queue.top();
		Queue.pop();
		CharacteristicInterval[1] = Queue.top();
		Queue.pop();
		CharacteristicInterval[2] = Queue.top();
		Queue.pop();
		CharacteristicInterval[3] = Queue.top();
		Queue.pop();

		#pragma omp parallel shared(CharacteristicInterval, testPoints) 
		{
			int numberOfThr = omp_get_thread_num();
			double tpk = 0.5 * (CharacteristicInterval[numberOfThr].pRight->x + CharacteristicInterval[numberOfThr].pLeft->x) - ((CharacteristicInterval[numberOfThr].pRight->z -
				CharacteristicInterval[numberOfThr].pLeft->z) / (2.0 * m));
			Points[numberOfThr].x = tpk;
			Points[numberOfThr].z = f(tpk);
		} 

		

		tpt = InsertUp_List(&testPoints, &Points[0]);
		Queue.push(Interval(ComputeR(*CharacteristicInterval[0].pLeft, *tpt, m), CharacteristicInterval[0].pLeft, tpt));
		Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval[0].pRight, m), tpt, CharacteristicInterval[0].pRight));

		tpt = InsertUp_List(&testPoints, &Points[1]);
		Queue.push(Interval(ComputeR(*CharacteristicInterval[1].pLeft, *tpt, m), CharacteristicInterval[1].pLeft, tpt));
		Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval[1].pRight, m), tpt, CharacteristicInterval[1].pRight));

		tpt = InsertUp_List(&testPoints, &Points[2]);
		Queue.push(Interval(ComputeR(*CharacteristicInterval[2].pLeft, *tpt, m), CharacteristicInterval[2].pLeft, tpt));
		Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval[2].pRight, m), tpt, CharacteristicInterval[2].pRight));

		tpt = InsertUp_List(&testPoints, &Points[3]);
		Queue.push(Interval(ComputeR(*CharacteristicInterval[3].pLeft, *tpt, m), CharacteristicInterval[3].pLeft, tpt));
		Queue.push(Interval(ComputeR(*tpt, *CharacteristicInterval[3].pRight, m), tpt, CharacteristicInterval[3].pRight));

		k += 4;


	} while (fabs((CharacteristicInterval[0].pRight->x - CharacteristicInterval[0].pLeft->x) > accuracy) || fabs((CharacteristicInterval[1].pRight->x - CharacteristicInterval[1].pLeft->x) > accuracy) 
		|| fabs((CharacteristicInterval[2].pRight->x - CharacteristicInterval[2].pLeft->x) > accuracy)
		|| fabs((CharacteristicInterval[3].pRight->x - CharacteristicInterval[3].pLeft->x) > accuracy) && ((k - 2) < iterations));

	itLeft = testPoints.begin();
	GlobalMin = itLeft->z;
	DotOfGlobalMin = itLeft->x;
	itLeft++;

	while (itLeft != testPoints.end())
	{
		if (GlobalMin > itLeft->z)
		{
			GlobalMin = itLeft->z;
			DotOfGlobalMin = itLeft->x;
		}

		itLeft++;
	}
	timeEnd = clock();
	
	cout.precision(9);
	cout << " Global minimum is: " << std::fixed << GlobalMin << endl;
	cout << " Dot of global minimum is: " << DotOfGlobalMin << endl;
	cout << " Number of experiments is: " << k << endl;
	cout << "Time of algorithm is: " << (timeEnd - timeStart) / CLOCKS_PER_SEC << endl;

}

