#ifndef BTC_H_
#define BTC_H_


struct Transport;

struct BTC{
	int		 n_BTC_elms;	  // Number of BTC elements
	int		*BTC_elm;	  // List of BTC elements
};


void output_transport_init_BTC(struct Transport *transport);
void output_transport_time_BTC(struct Transport *transport, double time);
void btc_check(struct Transport *transport);

#endif /* BTC_H_ */
