
histogram.o:     formato del fichero elf64-x86-64


Desensamblado de la sección .text:

0000000000000000 <histogram_input_init>:
  // choose the best version at runtime. If other compiler
  // is used, the default scalar version is used.
  h_input->priv_use_simd = 0;

#ifdef __GNUC__
  if (__builtin_cpu_supports("sse4.2") && __builtin_cpu_supports("popcnt")) {
   0:	8b 05 00 00 00 00    	mov    0x0(%rip),%eax        # 6 <histogram_input_init+0x6>
  // Check if SSE4.2 and POPCNT CPU extensions are supported.
  //
  // If the application is compiled using GCC, it will
  // choose the best version at runtime. If other compiler
  // is used, the default scalar version is used.
  h_input->priv_use_simd = 0;
   6:	c6 47 20 00          	movb   $0x0,0x20(%rdi)

#ifdef __GNUC__
  if (__builtin_cpu_supports("sse4.2") && __builtin_cpu_supports("popcnt")) {
   a:	f6 c4 01             	test   $0x1,%ah
   d:	74 08                	je     17 <histogram_input_init+0x17>
   f:	a8 04                	test   $0x4,%al
  11:	74 04                	je     17 <histogram_input_init+0x17>
    h_input->priv_use_simd = 1;
  13:	c6 47 20 01          	movb   $0x1,0x20(%rdi)
  }
#endif // __GNUC__

  h_input->in_read = (fastq_read_t*)read;
  17:	48 89 37             	mov    %rsi,(%rdi)
  1a:	c3                   	retq   
  1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000000020 <__histogram_count_cg>:
}

//------------------------------------------------------------------------------------

void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  20:	48 8b 07             	mov    (%rdi),%rax
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  23:	48 63 70 18          	movslq 0x18(%rax),%rsi
  27:	48 85 f6             	test   %rsi,%rsi
  2a:	74 49                	je     75 <__histogram_count_cg+0x55>
  2c:	4c 8b 40 08          	mov    0x8(%rax),%r8

//------------------------------------------------------------------------------------

void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;
  30:	31 c9                	xor    %ecx,%ecx

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  32:	31 c0                	xor    %eax,%eax

//------------------------------------------------------------------------------------

void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;
  34:	45 31 c9             	xor    %r9d,%r9d
  37:	eb 1c                	jmp    55 <__histogram_count_cg+0x35>
  39:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
		if (read->sequence[i] == 'C') {
			++nc;
		} else if (read->sequence[i] == 'G') {
			++ng;
  40:	80 fa 47             	cmp    $0x47,%dl
  43:	0f 94 c2             	sete   %dl
void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  46:	48 83 c0 01          	add    $0x1,%rax
		if (read->sequence[i] == 'C') {
			++nc;
		} else if (read->sequence[i] == 'G') {
			++ng;
  4a:	0f b6 d2             	movzbl %dl,%edx
  4d:	48 01 d1             	add    %rdx,%rcx
void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  50:	48 39 f0             	cmp    %rsi,%rax
  53:	73 17                	jae    6c <__histogram_count_cg+0x4c>
		if (read->sequence[i] == 'C') {
  55:	41 0f b6 14 00       	movzbl (%r8,%rax,1),%edx
  5a:	80 fa 43             	cmp    $0x43,%dl
  5d:	75 e1                	jne    40 <__histogram_count_cg+0x20>
void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  5f:	48 83 c0 01          	add    $0x1,%rax
		if (read->sequence[i] == 'C') {
			++nc;
  63:	49 83 c1 01          	add    $0x1,%r9
void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
  67:	48 39 f0             	cmp    %rsi,%rax
  6a:	72 e9                	jb     55 <__histogram_count_cg+0x35>
			++ng;
		}
	}

  // Store the occurrence count
  h_input->out_nc = nc;
  6c:	4c 89 4f 08          	mov    %r9,0x8(%rdi)
  h_input->out_ng = ng;
  70:	48 89 4f 10          	mov    %rcx,0x10(%rdi)
  74:	c3                   	retq   

//------------------------------------------------------------------------------------

void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;
  75:	31 c9                	xor    %ecx,%ecx
  77:	45 31 c9             	xor    %r9d,%r9d
  7a:	eb f0                	jmp    6c <__histogram_count_cg+0x4c>
  7c:	0f 1f 40 00          	nopl   0x0(%rax)

0000000000000080 <__histogram_count_cg_simd>:
  h_input->out_ng = ng;
}

//------------------------------------------------------------------------------------

void __histogram_count_cg_simd(histogram_input_t* h_input) {
  80:	55                   	push   %rbp
  81:	53                   	push   %rbx
  82:	48 83 ec 58          	sub    $0x58,%rsp
  const fastq_read_t* read = h_input->in_read;
  86:	48 8b 2f             	mov    (%rdi),%rbp
  h_input->out_ng = ng;
}

//------------------------------------------------------------------------------------

void __histogram_count_cg_simd(histogram_input_t* h_input) {
  89:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
  90:	00 00 
  92:	48 89 44 24 48       	mov    %rax,0x48(%rsp)
  97:	31 c0                	xor    %eax,%eax
  // be counted in multiples of 16 characters
  size_t vec_count = read->length - (read->length % 16);

  // Store the reference strings for comparing with all C's
  // and all G's in the SIMD registers
  const char ref_c[16] __attribute__ ((aligned (16))) = {
  99:	66 0f 6f 1d 00 00 00 	movdqa 0x0(%rip),%xmm3        # a1 <__histogram_count_cg_simd+0x21>
  a0:	00 

  // Calculate how many bases will be counted using vector
  // instructions and how many using a scalar loop, since
  // in the case of using vector instructions these must
  // be counted in multiples of 16 characters
  size_t vec_count = read->length - (read->length % 16);
  a1:	44 8b 5d 18          	mov    0x18(%rbp),%r11d
  const char ref_c[16] __attribute__ ((aligned (16))) = {
    'C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C'
  };

  const char ref_g[16] __attribute__ ((aligned (16))) = {
  a5:	66 0f 6f 15 00 00 00 	movdqa 0x0(%rip),%xmm2        # ad <__histogram_count_cg_simd+0x2d>
  ac:	00 

  // Calculate how many bases will be counted using vector
  // instructions and how many using a scalar loop, since
  // in the case of using vector instructions these must
  // be counted in multiples of 16 characters
  size_t vec_count = read->length - (read->length % 16);
  ad:	44 89 da             	mov    %r11d,%edx
  b0:	44 89 db             	mov    %r11d,%ebx
  b3:	c1 fa 1f             	sar    $0x1f,%edx
  b6:	c1 ea 1c             	shr    $0x1c,%edx
  b9:	41 8d 04 13          	lea    (%r11,%rdx,1),%eax

  // Store the reference strings for comparing with all C's
  // and all G's in the SIMD registers
  const char ref_c[16] __attribute__ ((aligned (16))) = {
  bd:	66 0f 7f 5c 24 20    	movdqa %xmm3,0x20(%rsp)

  // Calculate how many bases will be counted using vector
  // instructions and how many using a scalar loop, since
  // in the case of using vector instructions these must
  // be counted in multiples of 16 characters
  size_t vec_count = read->length - (read->length % 16);
  c3:	83 e0 0f             	and    $0xf,%eax
  c6:	29 d0                	sub    %edx,%eax
  c8:	29 c3                	sub    %eax,%ebx
  ca:	4c 63 cb             	movslq %ebx,%r9
  // Count the C and G in multiples of 16
  // Settings for the string comparation:
  // - Assume 8 bit unsigned integer characters.
  // - Do the comparation character by character.
  // - Return a mask with a bit per character.
  for (size_t i = 0; i < vec_count; i += 16) {
  cd:	4d 85 c9             	test   %r9,%r9
  const char ref_c[16] __attribute__ ((aligned (16))) = {
    'C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C'
  };

  const char ref_g[16] __attribute__ ((aligned (16))) = {
  d0:	66 0f 7f 54 24 30    	movdqa %xmm2,0x30(%rsp)
  // Count the C and G in multiples of 16
  // Settings for the string comparation:
  // - Assume 8 bit unsigned integer characters.
  // - Do the comparation character by character.
  // - Return a mask with a bit per character.
  for (size_t i = 0; i < vec_count; i += 16) {
  d6:	0f 84 c5 00 00 00    	je     1a1 <__histogram_count_cg_simd+0x121>
  dc:	31 c9                	xor    %ecx,%ecx

//------------------------------------------------------------------------------------

void __histogram_count_cg_simd(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;
  de:	31 f6                	xor    %esi,%esi
  e0:	45 31 c0             	xor    %r8d,%r8d
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpestrm (__m128i __X, int __LX, __m128i __Y, int __LY, const int __M)
{
  return (__m128i) __builtin_ia32_pcmpestrm128 ((__v16qi)__X, __LX,
  e3:	b8 10 00 00 00       	mov    $0x10,%eax
  e8:	4c 8b 55 08          	mov    0x8(%rbp),%r10
  ec:	0f 1f 40 00          	nopl   0x0(%rax)
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadu_si128 (__m128i const *__P)
{
  return (__m128i) __builtin_ia32_loaddqu ((char const *)__P);
  f0:	f3 41 0f 6f 0c 0a    	movdqu (%r10,%rcx,1),%xmm1
  f6:	89 c2                	mov    %eax,%edx
  f8:	66 0f 3a 60 d9 08    	pcmpestrm $0x8,%xmm1,%xmm3
  fe:	66 0f 7f 04 24       	movdqa %xmm0,(%rsp)
 103:	66 0f 3a 60 d1 08    	pcmpestrm $0x8,%xmm1,%xmm2

    // Perform the reduction on both results.
    // The reduction is performed on the lower 64 bits
    // of the SIMD register, which are unpacked with the
    // conversion intrinsic.
    nc += _mm_popcnt_u32(mask_nc[0]);
 109:	48 8b 14 24          	mov    (%rsp),%rdx
  // Count the C and G in multiples of 16
  // Settings for the string comparation:
  // - Assume 8 bit unsigned integer characters.
  // - Do the comparation character by character.
  // - Return a mask with a bit per character.
  for (size_t i = 0; i < vec_count; i += 16) {
 10d:	48 83 c1 10          	add    $0x10,%rcx

/* Calculate a number of bits set to 1.  */
extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_popcnt_u32 (unsigned int __X)
{
  return __builtin_popcount (__X);
 111:	f3 0f b8 d2          	popcnt %edx,%edx

    // Perform the reduction on both results.
    // The reduction is performed on the lower 64 bits
    // of the SIMD register, which are unpacked with the
    // conversion intrinsic.
    nc += _mm_popcnt_u32(mask_nc[0]);
 115:	48 63 d2             	movslq %edx,%rdx
 118:	49 01 d0             	add    %rdx,%r8
 11b:	66 0f 7f 44 24 10    	movdqa %xmm0,0x10(%rsp)
    ng += _mm_popcnt_u32(mask_ng[0]);
 121:	48 8b 54 24 10       	mov    0x10(%rsp),%rdx
 126:	f3 0f b8 d2          	popcnt %edx,%edx
 12a:	48 63 d2             	movslq %edx,%rdx
 12d:	48 01 d6             	add    %rdx,%rsi
  // Count the C and G in multiples of 16
  // Settings for the string comparation:
  // - Assume 8 bit unsigned integer characters.
  // - Do the comparation character by character.
  // - Return a mask with a bit per character.
  for (size_t i = 0; i < vec_count; i += 16) {
 130:	49 39 c9             	cmp    %rcx,%r9
 133:	77 bb                	ja     f0 <__histogram_count_cg_simd+0x70>
    ng += _mm_popcnt_u32(mask_ng[0]);
	}

  // Count the remaining Cs and Gs outside of the
  // 16 character boundary
  for (int i = vec_count; i < read->length; i++) {
 135:	44 39 db             	cmp    %r11d,%ebx
 138:	7d 48                	jge    182 <__histogram_count_cg_simd+0x102>
 13a:	48 8b 55 08          	mov    0x8(%rbp),%rdx
 13e:	f7 d3                	not    %ebx
 140:	4a 8d 4c 0a 01       	lea    0x1(%rdx,%r9,1),%rcx
 145:	4a 8d 04 0a          	lea    (%rdx,%r9,1),%rax
 149:	42 8d 14 1b          	lea    (%rbx,%r11,1),%edx
 14d:	48 01 d1             	add    %rdx,%rcx
 150:	eb 1b                	jmp    16d <__histogram_count_cg_simd+0xed>
 152:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
    if (read->sequence[i] == 'C') {
			++nc;
		} else if (read->sequence[i] == 'G') {
			++ng;
 158:	80 fa 47             	cmp    $0x47,%dl
 15b:	0f 94 c2             	sete   %dl
 15e:	48 83 c0 01          	add    $0x1,%rax
 162:	0f b6 d2             	movzbl %dl,%edx
 165:	48 01 d6             	add    %rdx,%rsi
    ng += _mm_popcnt_u32(mask_ng[0]);
	}

  // Count the remaining Cs and Gs outside of the
  // 16 character boundary
  for (int i = vec_count; i < read->length; i++) {
 168:	48 39 c8             	cmp    %rcx,%rax
 16b:	74 15                	je     182 <__histogram_count_cg_simd+0x102>
    if (read->sequence[i] == 'C') {
 16d:	0f b6 10             	movzbl (%rax),%edx
 170:	80 fa 43             	cmp    $0x43,%dl
 173:	75 e3                	jne    158 <__histogram_count_cg_simd+0xd8>
 175:	48 83 c0 01          	add    $0x1,%rax
			++nc;
 179:	49 83 c0 01          	add    $0x1,%r8
    ng += _mm_popcnt_u32(mask_ng[0]);
	}

  // Count the remaining Cs and Gs outside of the
  // 16 character boundary
  for (int i = vec_count; i < read->length; i++) {
 17d:	48 39 c8             	cmp    %rcx,%rax
 180:	75 eb                	jne    16d <__histogram_count_cg_simd+0xed>
  }

  // Store the occurrence count
  h_input->out_nc = nc;
  h_input->out_ng = ng;
}
 182:	48 8b 44 24 48       	mov    0x48(%rsp),%rax
 187:	64 48 33 04 25 28 00 	xor    %fs:0x28,%rax
 18e:	00 00 
			++ng;
		}
  }

  // Store the occurrence count
  h_input->out_nc = nc;
 190:	4c 89 47 08          	mov    %r8,0x8(%rdi)
  h_input->out_ng = ng;
 194:	48 89 77 10          	mov    %rsi,0x10(%rdi)
}
 198:	75 0e                	jne    1a8 <__histogram_count_cg_simd+0x128>
 19a:	48 83 c4 58          	add    $0x58,%rsp
 19e:	5b                   	pop    %rbx
 19f:	5d                   	pop    %rbp
 1a0:	c3                   	retq   

//------------------------------------------------------------------------------------

void __histogram_count_cg_simd(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;
 1a1:	31 f6                	xor    %esi,%esi
 1a3:	45 31 c0             	xor    %r8d,%r8d
 1a6:	eb 8d                	jmp    135 <__histogram_count_cg_simd+0xb5>
  }

  // Store the occurrence count
  h_input->out_nc = nc;
  h_input->out_ng = ng;
}
 1a8:	e8 00 00 00 00       	callq  1ad <__histogram_count_cg_simd+0x12d>
 1ad:	0f 1f 00             	nopl   (%rax)

00000000000001b0 <histogram_apply>:

//------------------------------------------------------------------------------------

void histogram_apply(histogram_input_t* h_input) {
 1b0:	53                   	push   %rbx
 1b1:	48 89 fb             	mov    %rdi,%rbx
 1b4:	48 83 ec 10          	sub    $0x10,%rsp
  // Initialize histogram values
  h_input->out_nc = 0;
 1b8:	48 c7 47 08 00 00 00 	movq   $0x0,0x8(%rdi)
 1bf:	00 
  h_input->out_ng = 0;
 1c0:	48 c7 47 10 00 00 00 	movq   $0x0,0x10(%rdi)
 1c7:	00 
  h_input->out_ncg = 0.0f;
 1c8:	c7 47 18 00 00 00 00 	movl   $0x0,0x18(%rdi)
  h_input->out_ngc = 0.0f;
 1cf:	c7 47 1c 00 00 00 00 	movl   $0x0,0x1c(%rdi)

  // Count ocurrences of C and G, using the scalar or
  // vector version if extensions are available.
  //if (h_input->priv_use_simd) {
    __histogram_count_cg_simd(h_input);
 1d6:	e8 00 00 00 00       	callq  1db <histogram_apply+0x2b>
  //} else {
  //  __histogram_count_cg(h_input);
  //}

  // Calculate the proportion of C and G
  h_input->out_ncg = (h_input->out_nc + h_input->out_ng == 0) ? (0.5) : 
 1db:	48 8b 43 08          	mov    0x8(%rbx),%rax
 1df:	48 89 c2             	mov    %rax,%rdx
 1e2:	48 03 53 10          	add    0x10(%rbx),%rdx
 1e6:	74 38                	je     220 <histogram_apply+0x70>
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
 1e8:	48 85 c0             	test   %rax,%rax
 1eb:	78 5b                	js     248 <histogram_apply+0x98>
 1ed:	48 85 d2             	test   %rdx,%rdx
 1f0:	f2 48 0f 2a c8       	cvtsi2sd %rax,%xmm1
 1f5:	78 6b                	js     262 <histogram_apply+0xb2>
 1f7:	f2 48 0f 2a c2       	cvtsi2sd %rdx,%xmm0
 1fc:	f2 0f 5e c8          	divsd  %xmm0,%xmm1
 200:	f3 0f 10 05 00 00 00 	movss  0x0(%rip),%xmm0        # 208 <histogram_apply+0x58>
 207:	00 
  //} else {
  //  __histogram_count_cg(h_input);
  //}

  // Calculate the proportion of C and G
  h_input->out_ncg = (h_input->out_nc + h_input->out_ng == 0) ? (0.5) : 
 208:	66 0f 14 c9          	unpcklpd %xmm1,%xmm1
 20c:	66 0f 5a c9          	cvtpd2ps %xmm1,%xmm1
 210:	f3 0f 5c c1          	subss  %xmm1,%xmm0
 214:	eb 15                	jmp    22b <histogram_apply+0x7b>
 216:	66 2e 0f 1f 84 00 00 	nopw   %cs:0x0(%rax,%rax,1)
 21d:	00 00 00 
 220:	f3 0f 10 05 00 00 00 	movss  0x0(%rip),%xmm0        # 228 <histogram_apply+0x78>
 227:	00 
 228:	0f 28 c8             	movaps %xmm0,%xmm1
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
  h_input->out_ngc = 1.0 - h_input->out_ncg;

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
 22b:	83 3d 00 00 00 00 01 	cmpl   $0x1,0x0(%rip)        # 232 <histogram_apply+0x82>
  //} else {
  //  __histogram_count_cg(h_input);
  //}

  // Calculate the proportion of C and G
  h_input->out_ncg = (h_input->out_nc + h_input->out_ng == 0) ? (0.5) : 
 232:	f3 0f 11 4b 18       	movss  %xmm1,0x18(%rbx)
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
  h_input->out_ngc = 1.0 - h_input->out_ncg;
 237:	f3 0f 11 43 1c       	movss  %xmm0,0x1c(%rbx)

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
 23c:	7e 42                	jle    280 <histogram_apply+0xd0>
  LOG_DEBUG_F("========= VALUES Ncg %f =========\n", h_input->out_ncg);
  LOG_DEBUG_F("========= VALUES Ngc %f =========\n", h_input->out_ngc);
}
 23e:	48 83 c4 10          	add    $0x10,%rsp
 242:	5b                   	pop    %rbx
 243:	c3                   	retq   
 244:	0f 1f 40 00          	nopl   0x0(%rax)
  //  __histogram_count_cg(h_input);
  //}

  // Calculate the proportion of C and G
  h_input->out_ncg = (h_input->out_nc + h_input->out_ng == 0) ? (0.5) : 
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
 248:	48 89 c1             	mov    %rax,%rcx
 24b:	83 e0 01             	and    $0x1,%eax
 24e:	48 d1 e9             	shr    %rcx
 251:	48 09 c1             	or     %rax,%rcx
 254:	48 85 d2             	test   %rdx,%rdx
 257:	f2 48 0f 2a c9       	cvtsi2sd %rcx,%xmm1
 25c:	f2 0f 58 c9          	addsd  %xmm1,%xmm1
 260:	79 95                	jns    1f7 <histogram_apply+0x47>
 262:	48 89 d0             	mov    %rdx,%rax
 265:	83 e2 01             	and    $0x1,%edx
 268:	48 d1 e8             	shr    %rax
 26b:	48 09 d0             	or     %rdx,%rax
 26e:	f2 48 0f 2a c0       	cvtsi2sd %rax,%xmm0
 273:	f2 0f 58 c0          	addsd  %xmm0,%xmm0
 277:	eb 83                	jmp    1fc <histogram_apply+0x4c>
 279:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
  h_input->out_ngc = 1.0 - h_input->out_ncg;

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
 280:	48 8b 03             	mov    (%rbx),%rax
 283:	41 b9 00 00 00 00    	mov    $0x0,%r9d
 289:	41 b8 00 00 00 00    	mov    $0x0,%r8d
 28f:	b9 85 00 00 00       	mov    $0x85,%ecx
 294:	ba 00 00 00 00       	mov    $0x0,%edx
 299:	be 00 00 00 00       	mov    $0x0,%esi
 29e:	bf 01 00 00 00       	mov    $0x1,%edi
 2a3:	48 8b 40 08          	mov    0x8(%rax),%rax
 2a7:	48 89 04 24          	mov    %rax,(%rsp)
 2ab:	31 c0                	xor    %eax,%eax
 2ad:	e8 00 00 00 00       	callq  2b2 <histogram_apply+0x102>
  LOG_DEBUG_F("========= VALUES Ncg %f =========\n", h_input->out_ncg);
 2b2:	83 3d 00 00 00 00 01 	cmpl   $0x1,0x0(%rip)        # 2b9 <histogram_apply+0x109>
 2b9:	7f 83                	jg     23e <histogram_apply+0x8e>
 2bb:	f3 0f 10 43 18       	movss  0x18(%rbx),%xmm0
 2c0:	41 b9 00 00 00 00    	mov    $0x0,%r9d
 2c6:	41 b8 00 00 00 00    	mov    $0x0,%r8d
 2cc:	b9 86 00 00 00       	mov    $0x86,%ecx
 2d1:	ba 00 00 00 00       	mov    $0x0,%edx
 2d6:	be 00 00 00 00       	mov    $0x0,%esi
 2db:	0f 5a c0             	cvtps2pd %xmm0,%xmm0
 2de:	bf 01 00 00 00       	mov    $0x1,%edi
 2e3:	b8 01 00 00 00       	mov    $0x1,%eax
 2e8:	e8 00 00 00 00       	callq  2ed <histogram_apply+0x13d>
  LOG_DEBUG_F("========= VALUES Ngc %f =========\n", h_input->out_ngc);
 2ed:	83 3d 00 00 00 00 01 	cmpl   $0x1,0x0(%rip)        # 2f4 <histogram_apply+0x144>
 2f4:	0f 8f 44 ff ff ff    	jg     23e <histogram_apply+0x8e>
 2fa:	f3 0f 10 43 1c       	movss  0x1c(%rbx),%xmm0
}
 2ff:	48 83 c4 10          	add    $0x10,%rsp
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
  h_input->out_ngc = 1.0 - h_input->out_ncg;

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
  LOG_DEBUG_F("========= VALUES Ncg %f =========\n", h_input->out_ncg);
  LOG_DEBUG_F("========= VALUES Ngc %f =========\n", h_input->out_ngc);
 303:	41 b9 00 00 00 00    	mov    $0x0,%r9d
}
 309:	5b                   	pop    %rbx
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
  h_input->out_ngc = 1.0 - h_input->out_ncg;

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
  LOG_DEBUG_F("========= VALUES Ncg %f =========\n", h_input->out_ncg);
  LOG_DEBUG_F("========= VALUES Ngc %f =========\n", h_input->out_ngc);
 30a:	0f 5a c0             	cvtps2pd %xmm0,%xmm0
 30d:	41 b8 00 00 00 00    	mov    $0x0,%r8d
 313:	b9 87 00 00 00       	mov    $0x87,%ecx
 318:	ba 00 00 00 00       	mov    $0x0,%edx
 31d:	be 00 00 00 00       	mov    $0x0,%esi
 322:	bf 01 00 00 00       	mov    $0x1,%edi
 327:	b8 01 00 00 00       	mov    $0x1,%eax
 32c:	e9 00 00 00 00       	jmpq   331 <histogram_apply+0x181>
