//(currently) true if length is a power of 2,3,5,7,11,13
inline bool IsASupportedLength( size_t length )
{
    while( length > 1 )
    {
        if( length % 2 == 0 )
            length /= 2;
        else if( length % 3 == 0 )
            length /= 3;
        else if( length % 5 == 0 )
            length /= 5;
        else if( length % 7 == 0 )
            length /= 7;
        else if (length % 11 == 0)
            length /= 11;
        else if (length % 13 == 0)
            length /= 13;
        else
            return false;
        }
    return true;
}

inline size_t FindBlue(size_t len)
{
    size_t p = 1;
    while(p < len)
        p <<= 1;
    return 2 * p;
}

enum oclfftType_
{
    OCLFFT_DIRECT_CLFFT = 1,      /* Lengths in all dimensions (1d, 2d, 3d) are directly supported by clFFT */
    OCLFFT_FULL_1D_CHIRP,         /* Lengths in all dimensions (1d) are unsupported by clFFT */
    OCLFFT_2D_CHIRP_X,            /* Lengths in x-dimension (2d) are unsupported by clFFT */
    OCLFFT_2D_CHIRP_Y,            /* Lengths in y-dimension (2d) are unsupported by clFFT */
    OCLFFT_FULL_2D_CHIRP,         /* Lengths in all dimensions (2d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_X,            /* Lengths in x-dimension (3d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_Y,            /* Lengths in y-dimension (3d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_Z,            /* Lengths in z-dimension (3d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_XY,           /* Lengths in x- and y-dimensions (3d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_XZ,           /* Lengths in x- and z-dimensions (3d) are unsupported by clFFT */
    OCLFFT_3D_CHIRP_YZ,           /* Lengths in y- and z-dimensions (3d) are unsupported by clFFT */
    OCLFFT_FULL_3D_CHIRP,         /* Lengths in all dimensions (3d) are unsupported by clFFT */
    ENDOCLFFTTYPES                /* The last value of the enum, and marks the length of oclfftType */
} oclfftType;
