#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ActivationFunction.h>
#include <matrix.h>
#include <ann.h>
#include <getopt.h>
#include <string.h>

int TRAINING_MODE = 1;
static char *ANNfile = "data/ann.dat";
static char *INPUTfile = "data/Input.dat";
static char *NN_OUTPUTfile = "data/NNOutput.dat";

struct option options[] = {
        { "usage",      0, 0, 'u'  },
        { "ann",    	1, 0, 'A'  },
        { "input",      1, 0, 'I'  },
        { "nn_output",  1, 0, 'N'  }
};
void print_usage( void ){

        fprintf( stdout, "\nUSAGE OPTIONS:\n"
                "  -u, --usage\n"
                "       print usage\n"
                "  -A, --ann\n"
                "       filename of initialized neural network (default data/ann.dat)\n"
                "  -I, --input\n"
                "       filename of network training input (default data/Input.dat)\n"
                "  -N, --nn_output\n"
                "       filename where to save trained neural network output (default data/NNOutput.dat)\n"
         );
        exit(0);

}

int main( int argc , char **argv ){

        ANN net;
	matrix INPUT, NN_OUTPUT;
        int i, j, N_sample;
        int opt;
        extern char *optarg;


        /* 0. Training options */
        do {
                opt = getopt_long( argc, argv, "uA:I:T:", options, NULL  );
                switch( opt ) {
                case 'u':       print_usage();
                                break;
                case 'A':       ANNfile = strdup( optarg );
                                break;
                case 'I':       INPUTfile = strdup( optarg );
                                break;
                case 'N':       NN_OUTPUTfile = strdup( optarg );
                                break;
                default:        break;
                }
        }while( opt >= 0  );

        /* Artificial Neural Network inizialization*/
        printf("LOADING DATA...\n");
        if( ANN_init( &net, ANNfile ) ){
		fprintf( stderr, "Initialization error\n" );
                return 1;
        }
        /* Input data acquisition*/
        N_sample = 0;
        if( ANN_DataRead( &INPUT, &N_sample, INPUTfile ) ){
		fprintf( stderr, "Data input acquisition error\n" );
                return 1;
        }

	
	if( matrix_init( &NN_OUTPUT, N_sample, net.N_output ) ){
		fprintf( stderr, "MAtrix initailization error\n" );
		return 1;
	}

        ANN_write( &net, stdout, W_A_TEXT);
	fprintf( stdout, "SIMULATION....\n" );
        for( i=0; i<N_sample; i++ ){
                /* aggiorno il vettore degli ingressi */
		for( j=0; j<net.N_input; j++ ){
                        net.input.vec[j] = INPUT.mat[i][j];
                }
		/* simulo la rete */
                if( ANN_sim( &net, &net.input, &net.output) ){
			fprintf( stderr, "Network simulation error\n" );
                        return 1;
                }
		/* aggiorno la matrice delle uscite */
                for( j=0; j<net.N_output; j++ ){
                        NN_OUTPUT.mat[i][j] = net.output.vec[j];
                }

        }

	ANN_DataWrite( &NN_OUTPUT, NN_OUTPUTfile );

        /* dynamic memory free*/
        ANN_destroy( &net );
	matrix_destroy( &INPUT );
	matrix_destroy( &NN_OUTPUT );
        
	return 0;
}

