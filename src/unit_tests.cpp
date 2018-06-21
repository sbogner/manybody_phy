#include "catch.hpp"
#include "system.h"

TEST_CASE("TEST: Pairing Model States"){

	CFermionSystem PairingModel(4, 4, 1.0, 0.5);

	REQUIRE(PairingModel.mb_states_[0][0] == 0);
	REQUIRE(PairingModel.mb_states_[0][1] == 1);

	REQUIRE(PairingModel.mb_states_[1][0] == 0);
	REQUIRE(PairingModel.mb_states_[1][1] == -1);

	REQUIRE(PairingModel.mb_states_[2][0] == 1);
	REQUIRE(PairingModel.mb_states_[2][1] == 1);

	REQUIRE(PairingModel.mb_states_[3][0] == 1);
	REQUIRE(PairingModel.mb_states_[3][1] == -1);

	REQUIRE(PairingModel.mb_states_[4][0] == 2);
	REQUIRE(PairingModel.mb_states_[4][1] == 1);

	REQUIRE(PairingModel.mb_states_[5][0] == 2);
	REQUIRE(PairingModel.mb_states_[5][1] == -1);

	REQUIRE(PairingModel.mb_states_[6][0] == 3);
	REQUIRE(PairingModel.mb_states_[6][1] == 1);

	REQUIRE(PairingModel.mb_states_[7][0] == 3);
	REQUIRE(PairingModel.mb_states_[7][1] == -1);
}