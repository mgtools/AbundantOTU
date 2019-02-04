///////////////////////////////////////////////////////////////
// this class is used to store a read
// by yongzhao
///////////////////////////////////////////////////////////////


#ifndef __READ_H__
#define __READ_H__

#include <vector>
#include <string>
#include <iostream>
#include <boost/foreach.hpp>
using namespace std;
using namespace boost;


// this class store the read and it's flow value
class CRead
{
public:
	vector<double> m_vdFlow;
	string m_sName;
	string m_sRead;

public:
	void PrintFlow()
	{
		BOOST_FOREACH(double d, m_vdFlow)
		{
			cout << d << " ";
		}
		cout << endl;
	}
};


#endif
