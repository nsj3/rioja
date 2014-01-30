#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

integer f_open(olist *);
integer s_rsle(cilist *);
integer do_fio(ftnint *, char *, ftnlen);
integer e_rsfe(void);
integer s_rsle(cilist *);
integer do_lio(ftnint *, ftnint *, char *, ftnlen);
integer e_rsle(void);
void s_copy(char *, char *, ftnlen, ftnlen);
integer f_clos(cllist *);
integer s_rsfe(cilist *);

int closef_(int chan)
{
   cllist cl;
   cl.cerr = 0;
   cl.cunit = chan;
   cl.csta = 0;
   return f_clos(&cl);
}
/* Subroutine */ int openf_(char *fname, char *title, char *format, int *ncoup, int *chan, int *tag, int *filenamelength, cllist *cl)
{
/*  ftnlen fname_len = *filenamelength;
    ftnlen title_len = 80;
   ftnlen format_len = 68;
*/   
    /* Format strings */
    static char fmt_1000[] = "(a80)";
    static char fmt_1001[] = "(a68,i2)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
/*    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_rsle(), do_lio(), 
	    e_rsle(), f_clos();
*/	    

    /* Subroutine */ 
/*     int s_copy(); */

    /* Local variables */
    extern integer len_trim__();
    static integer n;

    /* Fortran I/O blocks */
    static cilist io___1 = { 1, 0, 1, fmt_1000, 0 };
    static cilist io___2 = { 1, 0, 1, fmt_1001, 0 };
    static cilist io___3 = { 1, 0, 1, 0, 0 };



/*     open file, get title, format, ncoup */
/*     set tag to 0 if success, 1 if fail open */
/*                              2 if fail read */
/*                              3 eof */
/*                              4 ncoup = 0 */

    o__1.oerr = 1;
    o__1.ounit = *chan;
/*    o__1.ofnmlen = 80;     SJ Modified from f2c */
    o__1.ofnmlen = *filenamelength;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L12;
    }
    io___1.ciunit = *chan;
    i__1 = s_rsfe(&io___1);

    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_fio(&c__1, title, 80L);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsfe();
L100001:
    if (i__1 < 0) {
	goto L14;
    }
    if (i__1 > 0) {
	goto L13;
    }
    *ncoup = 0;
    io___2.ciunit = *chan;
    i__1 = s_rsfe(&io___2);
    if (i__1 != 0) {
	goto L100002;
    }
    i__1 = do_fio(&c__1, format, 68L);
    if (i__1 != 0) {
	goto L100002;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncoup), (ftnlen)sizeof(int));
    if (i__1 != 0) {
	goto L100002;
    }
    i__1 = e_rsfe();
L100002:
    if (i__1 < 0) {
	goto L14;
    }
    if (i__1 > 0) {
	goto L13;
    }
    if (*ncoup == 0) {
	io___3.ciunit = *chan;
	i__1 = s_rsle(&io___3);
	if (i__1 != 0) {
	    goto L100003;
	}
	i__1 = do_lio(&c__2, &c__1, (char *)&(*ncoup), (ftnlen)sizeof(int));
	if (i__1 != 0) {
	    goto L100003;
	}
	i__1 = e_rsle();
L100003:
	if (i__1 < 0) {
	    goto L14;
	}
	if (i__1 > 0) {
	    goto L13;
	}
    }
    if (*ncoup == 0) {
	goto L15;
    }

/*  add null terminator */

/*    n = len_trim__(title, 80L); */
    n = 79;
    i__1 = n;
    s_copy(title + i__1, "\000", n + 1 - i__1, 1L);
/*    n = len_trim__(format, 60L); */
    n=67;
    i__1 = n;
    s_copy(format + i__1, "\000", n + 1 - i__1, 1L);

    *tag = 0;
    
    return 0;
L12:
    *tag = 1;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L13:
    *tag = 2;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L14:
    *tag = 3;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L15:
    *tag = 4;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* openf_ */

/* Subroutine */ int getlin_(int *chan, char *format, int *ncoup, int *cursam, int *sp, double *abun, int *tag)
{

/* ftnlen format_len = *formatlen; */
    /* System generated locals */
    integer i__1, i__2;
    cilist ci__1;
    cllist cl__1;

    /* Builtin functions */
/*    integer s_rsfe(), do_fio(), e_rsfe(), f_clos(); */

    /* Local variables */
    static integer i;


/*      reads next current sample from cornell condensed file */
/*      tag = 0 = success, 1 = error in data */
/*                         2 = unexpected EOF */

    /* Parameter adjustments */
    --abun;
    --sp;

    /* Function Body */
    ci__1.cierr = 1;
    ci__1.ciend = 1;
    ci__1.ciunit = *chan;
    ci__1.cifmt = format;
    i__1 = s_rsfe(&ci__1);
    if (i__1 != 0) {
	goto L100004;
    }
    i__1 = do_fio(&c__1, (char *)&(*cursam), (ftnlen)sizeof(int));
    if (i__1 != 0) {
	goto L100004;
    }
    i__2 = *ncoup;
    for (i = 1; i <= i__2; ++i) {
	i__1 = do_fio(&c__1, (char *)&sp[i], (ftnlen)sizeof(int));
	if (i__1 != 0) {
	    goto L100004;
	}
	i__1 = do_fio(&c__1, (char *)&abun[i], (ftnlen)sizeof(double));
	if (i__1 != 0) {
	    goto L100004;
	}
    }
    i__1 = e_rsfe();
L100004:
    if (i__1 < 0) {
	goto L20;
    }
    if (i__1 > 0) {
	goto L21;
    }
    *tag = 0;
    return 0;
L20:
    *tag = 2;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L21:
    *tag = 1;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* getlin_ */

/* Subroutine */ int getl2_(int *chan, char *format, int *ncoup, int *cursam, double *abun, int *tag)
{
/* ftnlen format_len = *formatlen; */
    /* System generated locals */
    integer i__1, i__2;
    cilist ci__1;
    cllist cl__1;

    /* Builtin functions */
/*    integer s_rsfe(), do_fio(), e_rsfe(), f_clos(); */

    /* Local variables */
    static integer i;


/*      reads next current sample from cornell condensed file */
/*      tag = 0 = success, 1 = error in data */
/*                         2 = unexpected EOF */

    /* Parameter adjustments */
    --abun;

    /* Function Body */
    ci__1.cierr = 1;
    ci__1.ciend = 1;
    ci__1.ciunit = *chan;
    ci__1.cifmt = format;
    i__1 = s_rsfe(&ci__1);
    if (i__1 != 0) {
	goto L100005;
    }
    i__1 = do_fio(&c__1, (char *)&(*cursam), (ftnlen)sizeof(int));
    if (i__1 != 0) {
	goto L100005;
    }
    i__2 = *ncoup;
    for (i = 1; i <= i__2; ++i) {
	i__1 = do_fio(&c__1, (char *)&abun[i], (ftnlen)sizeof(double));
	if (i__1 != 0) {
	    goto L100005;
	}
    }
    i__1 = e_rsfe();
L100005:
    if (i__1 < 0) {
	goto L20;
    }
    if (i__1 > 0) {
	goto L21;
    }
    *tag = 0;
    return 0;
L20:
    *tag = 2;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L21:
    *tag = 1;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* getl2_ */

/* Subroutine */ int getnam_(int *chan, char *spnam, char *samnam, int *nsp, int *nsam, int *tag)
{
    /* Format strings */
    static char fmt_1002[] = "(80a1)";

    /* System generated locals */
    integer i__1, i__2;
    cllist cl__1;

    /* Builtin functions */
/*    integer s_rsfe(), do_fio(), e_rsfe(), f_clos(); */

    /* Local variables */
    static integer i;

    /* Fortran I/O blocks */
    static cilist io___7 = { 1, 0, 1, fmt_1002, 0 };
    static cilist io___9 = { 1, 0, 1, fmt_1002, 0 };



/*      reads next current sample from cornell condensed file */
/*      tag = 0 = success, 1 = unexpected end in species data */
/*                         2 = error in species data */
/*                         3 = unexpected end in sample data */
/*                         4 = error in sample data */

    /* Parameter adjustments */
    --samnam;
    --spnam;

    /* Function Body */
    io___7.ciunit = *chan;
    i__1 = s_rsfe(&io___7);
    if (i__1 != 0) {
	goto L100006;
    }
    i__2 = *nsp << 3;
    for (i = 1; i <= i__2; ++i) {
	i__1 = do_fio(&c__1, spnam + i, 1L);
	if (i__1 != 0) {
	    goto L100006;
	}
    }
    i__1 = e_rsfe();
L100006:
    if (i__1 < 0) {
	goto L20;
    }
    if (i__1 > 0) {
	goto L21;
    }
    io___9.ciunit = *chan;
    i__1 = s_rsfe(&io___9);
    if (i__1 != 0) {
	goto L100007;
    }
    i__2 = *nsam << 3;
    for (i = 1; i <= i__2; ++i) {
	i__1 = do_fio(&c__1, samnam + i, 1L);
	if (i__1 != 0) {
	    goto L100007;
	}
    }
    i__1 = e_rsfe();
L100007:
    if (i__1 < 0) {
	goto L22;
    }
    if (i__1 > 0) {
	goto L23;
    }
    *tag = 0;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L20:
    *tag = 1;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L21:
    *tag = 2;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L22:
    *tag = 3;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L23:
    *tag = 4;
    cl__1.cerr = 0;
    cl__1.cunit = *chan;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* getnam_ */

