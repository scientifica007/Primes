# مقارنة المشروع بمرجع خارجي احترافي: primesieve

تاريخ القياس: **31 أغسطس 2026**.

## الهدف

حتى هذه المرحلة كانت معظم benchmarks تقارن نسخ المشروع ببعضها. هذا مفيد هندسيًا لكنه لا يجيب:

> أين يقف التطبيق أمام implementation احترافي معروف؟

لذلك شغّلنا `primesieve` على **نفس GitHub Actions runner** وبـsingle thread، وفصلنا بين workloadين:

1. **count فقط**: أقصى throughput لـ`primesieve`، وليس مقارنة متكافئة مع كودنا.
2. **iteration + count + checksum + last prime**: C++ iterator يمر على كل عدد أولي، وهو أقرب إلى `consume_primes_compact_spec()` في مشروعنا.

التشغيل الرسمي:

```text
GitHub Actions run: 33346218939
commit: a29c0cdd5730ef5e76390905f32f2a868be638eb
N = 100,000,000
Python = 3.12.14
primesieve = 12.0
CPU = Intel Xeon 6973P-C
logical CPUs = 4
primesieve-reported L1 = 48 KiB
primesieve-reported L2 = 2048 KiB
L3 = 480 MiB
```

جميع الطرق أكدت:

```text
pi(100000000) = 5,761,455
last = 99,999,989
checksum = 279,209,790,387,276
```

---

## تطبيق المشروع

كل إعداد شُغّل مرتين، والرقم التالي هو الوسيط.

| التطبيق | segment | الزمن | algorithmic working set | strike attempts |
|---|---:|---:|---:|---:|
| Wheel-210 compact | 64 KiB | 13.266285 s | 71,335 B | 28,560,753 |
| Wheel-30030 compact | 1 MiB | **12.402345 s** | 1,125,697 B | 20,748,015 |

ملاحظة: هذه الأوقات أسرع من بعض تشغيلاتنا السابقة لأن GitHub hosted runners ليست عتادًا ثابتًا؛ لذلك لا تُقارن أزمنة من run مختلف كأنها على نفس CPU.

---

## primesieve: iteration workload

استخدمنا `primesieve::iterator` من C++، ومررنا على جميع الأوليات حتى `10^8` وحسبنا نفس:

```text
count + checksum + last
```

التشغيلان:

```text
0.008943868 s
0.008886538 s
```

الوسيط/المتوسط لنقطتين:

```text
~0.008915 s
```

مقارنة بهذا workload:

```text
Wheel-210 project / primesieve iterator   ~= 1488x
Wheel-30030 project / primesieve iterator ~= 1391x
```

أي أن تطبيق Python الحالي **ليس منافسًا أداءً** لـ`primesieve`.

لكن لا يجوز تفسير النسبة بأنها فرق خوارزمي صرف؛ فهي تجمع عدة عوامل:

- Python interpreter overhead؛
- loop per candidate / per strike في Python؛
- حسابات `divmod` وفهرسة عالية المستوى؛
- غياب SIMD؛
- غياب extreme loop unrolling؛
- غياب pre-sieving المتقدم؛
- غياب specialized small/medium/big prime paths؛
- غياب bucket scheduling؛
- اختلاف تمثيل wheel/sieve indexes؛
- C++ native compilation و`-O3 -march=native` في baseline.

إذن هذا benchmark يحدد **فجوة التطبيق**، ولا يثبت أن اختيار Wheel معين أبطأ رياضيًا بألف مرة.

---

## primesieve: count throughput

شغّلنا أيضًا:

```text
primesieve 100000000 --count --threads=1 --time
```

ثلاث مرات:

```text
0.007 s
0.006 s
0.007 s
```

والبرنامج اختار تلقائيًا:

```text
Sieve size = 384 KiB
Threads = 1
```

هذا يمثل سقف throughput لمهمة **العد فقط**، لذلك لا نقارنه مباشرة بدالة مشروعنا التي تمر على كل أولي وتحسب checksum.

---

## ماذا تعلمنا؟

### 1. التحسين داخل Python ليس benchmark بحثيًا كافيًا

أن تصبح Wheel-30030 أسرع من Wheel-210 في تطبيقنا لا يعني أننا تفوقنا على state of the art. أفضل نسخة لدينا ما زالت بعيدة جدًا عن تطبيق متخصص.

### 2. نتائج wheel/cache السابقة تبقى مفيدة، لكن بوصفها micro-engineering داخل مشروعنا

هي تجيب:

> ما التمثيل الأفضل ضمن هذا التصميم وهذا runtime؟

ولا تجيب:

> هل لدينا خوارزمية أسرع من المعروف؟

### 3. المقارنة الخوارزمية العادلة تحتاج port منخفض المستوى

إذا أردنا تقييم فكرة Wheel-30030 أو أي scheduling جديد على أنها **خوارزمية** لا **كود Python**، يجب تنفيذ hot loop في C/C++/Rust أو على الأقل extension native، ثم المقارنة single-thread مع `primesieve` على نفس workload.

### 4. Oliveira e Silva هو prior art لا baseline زمني مباشر في هذا القياس

خوارزمية bucket الخاصة بـOliveira e Silva صُممت خصوصًا للمجالات العالية عندما تصبح sieving primes الكبيرة مشكلة cache، و`primesieve` نفسه يستخدم هذا الأسلوب عند الحاجة للأعداد فوق 32 bit. عند `N=10^8 < 2^32` نحن لا نختبر السيناريو الذي صُممت له bucket scheduling بصورة مباشرة.

لذلك لا ننقل أرقام Oliveira e Silva التاريخية إلى جدولنا؛ كانت على عتاد ومجالات ومهمات مختلفة. نستخدمها كمرجع prior art، و`primesieve` كbaseline تنفيذي حديث.

---

## وضع المشروع بعد هذا القياس

### معروف/معاد تنفيذه

- segmented sieve؛
- wheel/primorial filtering؛
- bit packing؛
- cache-sized segments؛
- فكرة ضغط جداول wheel؛
- space/time tradeoffs العامة.

### غير منفذ لدينا بعد لكنه معروف

- bucket scheduling للأوليات الكبيرة؛
- specialized algorithms حسب حجم sieving prime؛
- SIMD pre-sieving؛
- branch-optimized/unrolled inner loops.

### ما زال يستحق البحث

1. الدراسة الرقمية عبر قواعد العد **بعد conditioning modular صارم** وnull models سليمة.
2. أي statistic يبقى ثابتًا عبر نطاقات مستقلة ولا تفسره النتائج المعروفة.
3. نموذج cost/cache يتنبأ باختيار wheel/segment عبر معماريات متعددة، إذا أثبت ميزة على heuristics المعروفة.
4. أي فكرة خوارزمية جديدة بعد تنفيذ prior art الأساسي أولًا ومقارنتها في لغة منخفضة المستوى.

---

## الملفات المرتبطة

- `PRIOR_ART.md` — خريطة الأعمال السابقة وتصنيف الجدة.
- `baseline_primesieve.cpp` — workload C++ المكافئ تقريبًا.
- `benchmark_project_external_baseline.py` — إعدادات المشروع المستخدمة.
- `.github/workflows/external-baseline.yml` — benchmark يدوي قابل لإعادة التشغيل.

الـartifact الأصلي للتشغيل `33346218939` احتوى:

```text
project-baseline.txt
primesieve-iterator-baseline.txt
primesieve-count-baseline.txt
external-baseline-environment.txt
```
