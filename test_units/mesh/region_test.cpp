/*
 * material_dispatch_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include <gtest/gtest.h>
#include "mesh/region.hh"

TEST(Region, all) {
    EXPECT_THROW( Region::db().size(), RegionDB::ExcSizeWhileOpen );

    {
    Region r=Region::db().add_region(1001,"top", true);
    EXPECT_EQ(0, r.idx() );
    EXPECT_TRUE( r.is_boundary() );
    EXPECT_EQ(0, r.boundary_idx() );
    EXPECT_EQ(0, Region::db().add_region(1001,"top", true).idx() );
    }

    {
    Region r=Region::db().add_region(1002,"inside 1", false);
    EXPECT_EQ(1, r.idx() );
    EXPECT_EQ(0, r.bulk_idx() );
    EXPECT_EQ("inside 1", r.label() );
    EXPECT_EQ(1002, r.id() );
    }

    {
    Region r=Region::db().add_region(1003,"inside 2", false);
    EXPECT_EQ(3, r.idx() );
    EXPECT_EQ(1, r.bulk_idx() );
    }

    {
    Region r=Region::db().add_region(1004,"bottom", true);
    EXPECT_EQ(2, r.idx() );
    EXPECT_EQ(1, r.boundary_idx() );
    EXPECT_EQ("bottom", r.label() );
    EXPECT_EQ(1004, r.id() );
    }

    Region::db().add_region(1005,"side", true);
    EXPECT_THROW( Region::db().add_region(1005,"new", false) , RegionDB::ExcInconsistentAdd);
    EXPECT_THROW( Region::db().add_region(1001,"bottom", false) , RegionDB::ExcInconsistentAdd);

    Region::db().close();

    EXPECT_EQ(6, Region::db().size());
    EXPECT_EQ(3, Region::db().boundary_size());
    EXPECT_EQ(2, Region::db().bulk_size());

    EXPECT_THROW( Region::db().add_region(1006,"side_", true) , RegionDB::ExcAddingIntoClosed );

}


