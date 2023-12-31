
typedef union  {
   Vec vec;
   VList *vecl;
   Flt flt;
   Flt *fltptr;
   Object *obj;
   NODE_PTR exper;
   LIST_PTR elist;
   char *name;
   void *data;
   csgnodeptr csgtree;
   Special_Surface *surf;
   Texture *text;
   texture_map_entries text_map;
   texture_fn_entries text_fn;
   map_entries cmap_entry;
   Transform *trns;
   ostackptr objlist;
   tstackptr textlist;
   Particle *part;
   Light *lgt;
} YYSTYPE;
extern YYSTYPE yylval;
#define ACCELERATION 257
#define ACOS 258
#define AMBIENT 259
#define AND_EXPER 260
#define ANGLE 261
#define ANTIALIAS 262
#define ANTIALIAS_THRESHOLD 263
#define APERTURE 264
#define ARRAY 265
#define ASSIGNMENT 266
#define ASIN 267
#define ASPECT 268
#define AT 269
#define ATAN 270
#define ATAN_TWO 271
#define AVOID 272
#define BACKGROUND 273
#define BEZIER 274
#define BIAS 275
#define BIRTH 276
#define BLINN 277
#define BLOB 278
#define BOUNDING_BOX 279
#define BOX 280
#define BRILLIANCE 281
#define BUMP_SCALE 282
#define CEIL 283
#define CHECKER 284
#define CHEIGHT_FIELD 285
#define CHEIGHT_FN 286
#define COLOR 287
#define COLOR_MAP 288
#define COLOR_WHEEL 289
#define CONCAT 290
#define CONCAVE 291
#define CONDITIONAL_EXPER 292
#define CONE 293
#define CONTOUR 294
#define COOK 295
#define COS 296
#define COSH 297
#define COUNT 298
#define CSG 299
#define CYLINDER 300
#define CYLINDRICAL_BUMPMAP 301
#define CYLINDRICAL_IMAGEMAP 302
#define CYLINDRICAL_INDEXED 303
#define DEATH 304
#define DEFINE 305
#define DEGREES 306
#define DEPTH 307
#define DEPTHMAPPED_LIGHT 308
#define DIFFUSE 309
#define DIRECTIONAL_LIGHT 310
#define DISC 311
#define DISPLACE 312
#define DITHER 313
#define DIV_EXPER 314
#define DNOISE 315
#define DOT_EXPER 316
#define DRAW 317
#define ELSE 318
#define END_FRAME 319
#define ENVIRONMENT 320
#define ENVIRONMENT_MAP 321
#define EQUAL_EXPER 322
#define EXP 323
#define EXPRESSION_SYM 324
#define FABS 325
#define FBM 326
#define FERRARI 327
#define FILE_FLUSH 328
#define FLARE 329
#define FLOCK 330
#define FLOOR 331
#define FNOISE 332
#define FOCAL_DISTANCE 333
#define FMOD 334
#define FRAME 335
#define FRAME_TIME 336
#define FREQUENCY 337
#define FROM 338
#define FUNCTION 339
#define GAIN 340
#define GAUSSIAN 341
#define GLYPH 342
#define GREATER_EXPER 343
#define GRIDDED 344
#define GTEQ_EXPER 345
#define HAZE 346
#define HEIGHT_FIELD 347
#define HEIGHT_FN 348
#define HEIGHT_MAP 349
#define HEXAGON 350
#define HITHER 351
#define HYPERTEXTURE 352
#define I_EXPER 353
#define IF 354
#define IMAGE 355
#define IMAGE_FORMAT 356
#define IMAGE_WINDOW 357
#define INCLUDE 358
#define INDEXED 359
#define INDEXED_MAP 360
#define LATHE 361
#define LAYERED 362
#define LEGENDRE 363
#define LENSES 364
#define LESS_EXPER 365
#define LIGHT 366
#define LN 367
#define LOG 368
#define LOOKUP_FUNCTION 369
#define LTEQ_EXPER 370
#define MAXT 371
#define MAX_SAMPLES 372
#define MAX_TRACE_DEPTH 373
#define MICROFACET 374
#define MINT 375
#define MINUS_EXPER 376
#define N_EXPER 377
#define NO_SHADOW 378
#define NOEVAL 379
#define NOISE 380
#define NORMAL 381
#define NOT_EXPER 382
#define NUM 383
#define NURB 384
#define OBJECT 385
#define OBJECT_SYM 386
#define OCTAVES 387
#define OR_EXPER 388
#define OPACITY 389
#define OUTFILE 390
#define P_EXPER 391
#define PARABOLA 392
#define PARAMETRIC 393
#define PARTICLE 394
#define PARTICLE_SYM 395
#define PATCH 396
#define PHASE 397
#define PHONG 398
#define PIXEL_ENCODING 399
#define PIXELSIZE 400
#define PLANAR_BUMPMAP 401
#define PLANAR_IMAGEMAP 402
#define PLANE 403
#define PLUS_EXPER 404
#define POINT 405
#define POLYGON 406
#define POLYNOMIAL 407
#define POSITION 408
#define POSITION_FUNCTION 409
#define POSITION_SCALE 410
#define POWER_EXPER 411
#define RADIANS 412
#define RAMP 413
#define RANDOM 414
#define RAW 415
#define REFLECT 416
#define REFLECTION 417
#define REITZ 418
#define RESOLUTION 419
#define RIPPLE 420
#define ROOT_SOLVER 421
#define ROTATE 422
#define SAWTOOTH 423
#define SCALE 424
#define SEED 425
#define SHADING_FLAGS 426
#define SHEAR 427
#define SIN 428
#define SINH 429
#define SIZE 430
#define SMOOTH_HEIGHT_FIELD 431
#define SMOOTH_HEIGHT_FN 432
#define SMOOTH_CHEIGHT_FIELD 433
#define SMOOTH_CHEIGHT_FN 434
#define SMOOTH_SHEIGHT_FIELD 435
#define SMOOTH_SHEIGHT_FN 436
#define SPACING 437
#define SPECIAL 438
#define SPECIAL_SURFACE_SYM 439
#define SHEIGHT_FIELD 440
#define SHEIGHT_FN 441
#define SPHERICAL_BUMPMAP 442
#define SPHERICAL_IMAGEMAP 443
#define SPHERICAL_INDEXED 444
#define SPLINE 445
#define SPOT_LIGHT 446
#define SQRT 447
#define SPECULAR 448
#define SPHERE 449
#define START_FRAME 450
#define STATIC 451
#define STRING 452
#define STURM 453
#define SUBSCRIPT_EXPER 454
#define SUMMED 455
#define SUPERQ 456
#define SURFACE 457
#define SURFACE_SYM 458
#define SYSTEM 459
#define SWEEP 460
#define TAN 461
#define TANH 462
#define TERM 463
#define TEXTURE 464
#define TEXTURE_MAP 465
#define TEXTURE_MAP_SYM 466
#define TEXTURE_SYM 467
#define TEXTURED_LIGHT 468
#define TIMES_EXPER 469
#define TOKEN 470
#define TORUS 471
#define TOTAL_FRAMES 472
#define TRACE 473
#define TRANSFORM 474
#define TRANSFORM_SYM 475
#define TRANSLATE 476
#define TRANSMISSION 477
#define TURBULENCE 478
#define UMINUS_EXPER 479
#define UP 480
#define U_EXPER 481
#define UU_EXPER 482
#define UV_EXPER 483
#define UW_EXPER 484
#define U_STEPS 485
#define V_STEPS 486
#define UV 487
#define UV_STEPS 488
#define UV_BOUNDS 489
#define VELOCITY 490
#define VIETA 491
#define VIEWPOINT 492
#define VISIBLE 493
#define VAL_EXPER 494
#define VEC_EXPER 495
#define VECTOR_EXPER 496
#define WAVE 497
#define W_EXPER 498
#define W_STEPS 499
#define X_EXPER 500
#define X_OFFSET 501
#define Y_EXPER 502
#define Y_OFFSET 503
#define YON 504
#define Z_EXPER 505
#define AND_SYM 506
#define OR_SYM 507
#define LTEQ_SYM 508
#define GTEQ_SYM 509
#define EQUAL_SYM 510
#define NEQUAL_SYM 511
#define UMINUS 512
