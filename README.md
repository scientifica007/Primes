# Primes

[![Tests](https://github.com/scientifica007/Primes/actions/workflows/tests.yml/badge.svg)](https://github.com/scientifica007/Primes/actions/workflows/tests.yml)

مشروع تجريبي لدراسة الأعداد الأولية من عدة زوايا: بناء الغربلة خطوةً بخطوة، تمثيل المرشحين، الذاكرة والـcache، عجلات primorial، وتغيير قواعد العد.

## الحالة البحثية — ابدأ من هنا

هذا المستودع **ليس حاليًا ادعاءً بخوارزمية جديدة للأعداد الأولية**. بعد مراجعة الأعمال السابقة، تبيّن أن معظم المسار الخوارزمي الذي وصلنا إليه — segmented sieve، wheel/primorial filtering، bit packing، cache-sized segments — له prior art معروف.

قبل متابعة أي تجربة جديدة اقرأ:

- **[RESEARCH_STATUS.md](RESEARCH_STATUS.md)** — الحالة الحالية وNovelty Gate.
- **[PRIOR_ART.md](PRIOR_ART.md)** — خريطة الأعمال السابقة: Singleton، Bays–Hudson، Pritchard، Sorenson، Oliveira e Silva، Járai–Vatai، Atkin–Bernstein و`primesieve`، إضافة إلى أدبيات الأنماط الرقمية.
- **[RESULTS_EXTERNAL_BASELINE.md](RESULTS_EXTERNAL_BASELINE.md)** — مقارنة مباشرة مع `primesieve` على نفس GitHub runner.
- **[RESULTS_BASE_REPRESENTATIONS.md](RESULTS_BASE_REPRESENTATIONS.md)** — دراسة قواعد العد حتى `N = 10^8`.
- **[RESULTS_COMPACT_CACHE.md](RESULTS_COMPACT_CACHE.md)** — دراسة wheel/cache والجداول المضغوطة.

### معيار الجدة من الآن فصاعدًا

لا يكفي أن تصبح نسخة أسرع من نسخة أقدم في هذا المستودع. قبل اعتبار نتيجة تقدمًا بحثيًا نسأل:

1. هل الفكرة موجودة في الأدبيات أو في تطبيقات احترافية مثل `primesieve`؟
2. هل هي مجرد إعادة صياغة لـmodular arithmetic؟
3. إذا كانت إحصائية، هل يوجد null model وتصحيح multiple testing وvalidation مستقل؟
4. إذا كانت ادعاء أداء، هل قورنت بbaseline خارجي، وهل فُصل أثر اللغة عن أثر الخوارزمية؟

## خط أساس خارجي

عند `N = 100,000,000`، على **نفس GitHub Actions runner** في التشغيل الرسمي `33346218939`:

```text
Project Wheel-210, Python, enumerate+checksum:   ~13.27 s
Project Wheel-30030, Python, enumerate+checksum: ~12.40 s
primesieve C++ iterator, enumerate+checksum:      ~0.0089 s
primesieve single-thread count only:              ~0.006–0.007 s
```

الفارق لا يُفسَّر كفرق خوارزمي صرف؛ تطبيقنا يعمل في Python ويخلو من عدد من تقنيات التنفيذ المنخفض المستوى الموجودة في `primesieve` مثل pre-sieving المتقدم، SIMD، loop unrolling، specialized prime paths وbucket scheduling. الغرض من هذا baseline هو منع المبالغة في تفسير تحسينات داخلية للمشروع.

## الفكرة الأصلية

بدأ المشروع من بناء غربلة بصورة قائمة/طرح:

1. نولد الأعداد الطبيعية حتى `N`.
2. نحذف الزوجيات.
3. نعالج `3` ثم `5` ثم بقية الأوليات.
4. بعد معالجة `p` نستفيد من حقيقة أن أي مركب له عامل أولي لا يتجاوز جذره.
5. تطورت الصياغة إلى إبقاء العدد الأولي نفسه وحذف مركباته فقط، والبدء من `p²`.

> العدد `1` ليس أوليًا ولا مركبًا.

هذه الصياغة التعليمية قادتنا تدريجيًا إلى بنى معروفة من غربال إراتوستينس، ثم إلى segmented/wheel implementations.

## التطبيقات الأساسية

### `original_method.py`

نسخة مرجعية تحافظ قدر الإمكان على بنية التجربة الأصلية: قوائم، مضاعفات، طرح، وتصنيف بين مربعات الأوليات.

```bash
python original_method.py 100
```

### `retain_prime_method.py`

صياغة مباشرة لفكرة **إبقاء العدد الأولي** داخل مجموعة العمل وحذف مضاعفاته المركبة، باستخدام `set`.

```bash
python retain_prime_method.py 100
```

التقرير: **[RESULTS_RETAIN_PRIME.md](RESULTS_RETAIN_PRIME.md)**

### `retain_prime_compact_method.py`

ضغط فضاء المرشحين:

- العدد `2` منفرد؛
- الأعداد الفردية فقط؛
- بت واحد لكل مرشح؛
- الشطب من `p²` وبخطوة `2p`.

عند `N = 100000` تحتاج مصفوفة المرشحين نفسها إلى `6250` بايت فقط.

التقرير: **[RESULTS_COMPACT_MEMORY.md](RESULTS_COMPACT_MEMORY.md)**

### `retain_prime_packed_output_method.py`

يضغط فضاء المرشحين والناتج معًا:

- bitset للمرشحين؛
- `array('I')`/`array('Q')` للناتج الكامل؛
- streaming عند عدم الحاجة للاحتفاظ بكل النتائج.

```bash
python retain_prime_packed_output_method.py 100
python retain_prime_packed_output_method.py 100 --stream
```

التقرير: **[RESULTS_PACKED_OUTPUT.md](RESULTS_PACKED_OUTPUT.md)**

### `segmented_method.py`

غربلة مجزأة: بدل bitset لكل المجال، يعالج مقطعًا ثابت السعة ثم يتخلص منه.

السعة الافتراضية الأصلية:

```text
32768 عددًا فرديًا = 4096 bytes للمقطع
```

التقرير: **[RESULTS_SEGMENTED.md](RESULTS_SEGMENTED.md)**

> ملاحظة prior art: segmented sieve معروف منذ أعمال مبكرة منها Singleton (1969) وBays–Hudson (1977). انظر `PRIOR_ART.md`.

### Wheel / Primorial sieves

طُورت تجارب:

```text
6 -> 30 -> 210 -> 2310 -> 30030
```

حيث لا نمثل إلا البواقي المتباينة أوليًا مع modulus العجلة.

ملفات وتقارير مرتبطة تشمل:

- `primorial_wheel_segmented_method.py`
- `compact_primorial_wheel_method.py`
- **[RESULTS_WHEEL30.md](RESULTS_WHEEL30.md)**
- **[RESULTS_WHEEL210.md](RESULTS_WHEEL210.md)**
- **[RESULTS_COMPACT_CACHE.md](RESULTS_COMPACT_CACHE.md)**

> ملاحظة prior art: wheel sieves وprimorial filtering معروفان خصوصًا من أعمال Paul Pritchard في 1981–1983.

### `compact_primorial_wheel_method.py`

نسخة cache-aware تضغط جداول العجلة إلى fixed-width arrays وتعيد استخدام `WheelSpec`.

في القياسات الرسمية انخفض payload جدول Wheel-30030 إلى نحو `71.6 KB` بدل مئات الكيلوبايت من كائنات Python في التمثيل القديم.

## دراسة قواعد العد

`base_representation_study.py` والملفات المرتبطة تدرس القواعد:

```text
b = 2, 3, ..., 100
```

مع اهتمام بـ:

```text
6, 30, 210, 2310
```

وتم التوسع حتى:

```text
N = 100,000,000
```

شملت الدراسة suffixes، entropy، transitions، digit sums، علاقات `b-1` و`b+1`، ومقارنة prime/composite مع null baselines.

النتيجة الحالية: معظم الإشارات القوية رجعت إلى modular structure معروف، ولم يظهر حتى الآن نمط رقمي قوي مستقل عنها.

التقرير: **[RESULTS_BASE_REPRESENTATIONS.md](RESULTS_BASE_REPRESENTATIONS.md)**

## Prior Art والاتجاهات المفتوحة

**[PRIOR_ART.md](PRIOR_ART.md)** يصنف كل اتجاه إلى:

- `KNOWN`
- `REPLICATION/ENGINEERING`
- `PROJECT-SPECIFIC BENCHMARK`
- `OPEN/EXPERIMENTAL`

أهم الاتجاهات التي ما زالت تستحق الفحص دون ادعاء جدة مسبق:

1. statistics رقمية بعد conditioning modular صارم، null models، multiple-testing correction وvalidation مستقل؛
2. نموذج cost/cache تنبؤي لاختيار wheel/segment عبر معماريات متعددة، إذا تفوق على heuristics المعروفة؛
3. أي فكرة غربلة جديدة فقط بعد تنفيذ/مقارنة prior art مثل bucket scheduling؛
4. port منخفض المستوى إذا كان الهدف مقارنة الخوارزمية نفسها بدل قياس overhead لغة Python.

## Baseline الاحترافي

الملفات:

- `baseline_primesieve.cpp`
- `benchmark_project_external_baseline.py`
- `.github/workflows/external-baseline.yml`

تسمح بإعادة مقارنة المشروع مع `primesieve` على نفس runner. الـworkflow يدوي حتى لا يعاد benchmark ثقيل مع كل commit.

التقرير: **[RESULTS_EXTERNAL_BASELINE.md](RESULTS_EXTERNAL_BASELINE.md)**

## التحقق

### `ci_verify.py`

مرجع تحقق مستقل يعتمد القسمة التجريبية حتى `sqrt(n)` ويقارن التطبيقات ضمن نطاق CI.

```bash
python ci_verify.py 4999
```

### اختبارات محلية

```bash
python -m unittest -v
```

تشمل الحالات الحدية، تطابق الطرق، bitsets، arrays، streaming، segmentation، wheels والجداول المضغوطة.

## GitHub Actions

### Tests

`.github/workflows/tests.yml` يعمل عند كل `push` إلى `main` وكل Pull Request على:

```text
Python 3.11
Python 3.12
Python 3.13
```

ويشغل compile، unit tests، والتحقق المستقل حتى `N = 4999`.

### Benchmarks

المسارات الثقيلة أصبحت يدوية في الحالات التي لا يلزم تشغيلها مع كل commit، ومنها:

- المقارنة العامة؛
- segmented benchmarks؛
- compact wheel/cache sweeps؛
- external `primesieve` baseline.

الـartifacts تحفظ بيئة التشغيل والنتائج بحيث لا نعتمد على زمن بلا CPU/commit موثق.

## النتائج الموثقة

- **[RESULTS.md](RESULTS.md)** — خط الأساس الأولي.
- **[RESULTS_RETAIN_PRIME.md](RESULTS_RETAIN_PRIME.md)** — إبقاء العدد الأولي.
- **[RESULTS_GITHUB_ACTIONS.md](RESULTS_GITHUB_ACTIONS.md)** — أول benchmark رسمي.
- **[RESULTS_COMPACT_MEMORY.md](RESULTS_COMPACT_MEMORY.md)** — ضغط فضاء المرشحين.
- **[RESULTS_PACKED_OUTPUT.md](RESULTS_PACKED_OUTPUT.md)** — ضغط الناتج والـstreaming.
- **[RESULTS_SEGMENTED.md](RESULTS_SEGMENTED.md)** — التجزئة والتوسع.
- **[RESULTS_WHEEL30.md](RESULTS_WHEEL30.md)** — Wheel-30.
- **[RESULTS_WHEEL210.md](RESULTS_WHEEL210.md)** — Wheel-210.
- **[RESULTS_BASE_REPRESENTATIONS.md](RESULTS_BASE_REPRESENTATIONS.md)** — قواعد العد حتى `10^8`.
- **[RESULTS_COMPACT_CACHE.md](RESULTS_COMPACT_CACHE.md)** — الجداول المضغوطة وحجم المقطع/cache.
- **[RESULTS_EXTERNAL_BASELINE.md](RESULTS_EXTERNAL_BASELINE.md)** — المقارنة مع `primesieve`.

## هدف المشروع بعد مراجعة prior art

الهدف لم يعد مجرد جعل نسخة المشروع أسرع من النسخة السابقة. الهدف الآن هو:

> **تمييز ما هو معروف وما هو إعادة تنفيذ وما هو سؤال مفتوح، ثم تصميم تجارب يمكن أن تفشل بوضوح وتنتج أدلة قابلة لإعادة الإنتاج.**

أي نتيجة مستقبلية يجب أن تمر عبر **[RESEARCH_STATUS.md](RESEARCH_STATUS.md)** و**[PRIOR_ART.md](PRIOR_ART.md)** قبل وصفها بأنها تقدم بحثي.
