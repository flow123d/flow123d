#ifndef REACTION_H_
#define REACTION_H_

//=============================================================================
// STRUCTURE OF THE REACTION
/*! @brief STRUCTURE OF THE REACTION (ONE SPECIE)
 *
 *
 */
//=============================================================================
struct Reaction
{
    int                     id; // reaction ID
    int                     sbi; // substance ID
    int                     type; // type of reaction
    double                  *coef; // type dependent coefficent set
};

void read_reaction_list( struct Transport *transport );
void transport_reaction(struct Transport *transport, int elm_pos, int sbi);

#endif /* REACTION_H_ */
