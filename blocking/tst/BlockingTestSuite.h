/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, write)
{
	gmds::blocking::Blocking bl;
	bl.createGrid();

	bl.writeMokaFile("aaa.moka");
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/