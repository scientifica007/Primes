# Primes

[![Tests](https://github.com/scientifica007/Primes/actions/workflows/tests.yml/badge.svg)](https://github.com/scientifica007/Primes/actions/workflows/tests.yml)

تجربة منهجية للبحث عن الأعداد الأولية ومقارنة صيغ مختلفة لتنفيذ الفكرة نفسها.

## الفكرة الأصلية

نبدأ بقائمة الأعداد الطبيعية حتى حد أعلى `N`، ثم:

1. نحذف الأعداد الزوجية، وتسمى القائمة الناتجة **أ**.
2. نضع `2` و`3` في قائمة الأعداد الأولية.
3. نحذف مضاعفات `3` من **أ**، وتسمى القائمة الناتجة **ب**.
4. كل عدد باقٍ يحقق `3 < n < 3²` هو عدد أولي.
5. نأخذ العدد الأولي التالي `5`، ونحذف مضاعفاته من **ب** لنحصل على قائمة جديدة.
6. كل عدد باقٍ يحقق `5 < n < 5²` هو أولي، مع عدم إعادة إدراج الأعداد التي صُنفت سابقًا.
7. ننتقل إلى العدد الأولي التالي ونكرر العملية.

الأساس الرياضي هو أن العدد المركب `n` يملك عاملًا أوليًا لا يتجاوز `sqrt(n)`. لذلك، بعد حذف مضاعفات جميع الأعداد الأولية حتى `p`، لا يمكن أن يبقى عدد مركب أصغر من `p²` لم يُحذف بعامل أصغر.

> العدد `1` مستبعد من التصنيف الأولي صراحةً لأنه ليس أوليًا ولا مركبًا.

## الملفات

### `original_method.py`

نسخة مرجعية تحافظ قدر الإمكان على بنية التجربة الأصلية: تنشئ القوائم، تنشئ المضاعفات، تطرحها، وتصنف الأعداد المتبقية بين مربعات الأعداد الأولية.

```bash
python original_method.py 100
```

### `retain_prime_method.py`

تجربة **إبقاء العدد الأولي** داخل مجموعة العمل وحذف مضاعفاته المركبة فقط. تستخدم `set` لتمثيل مباشر للفكرة.

```bash
python retain_prime_method.py 100
```

التقرير:

**[RESULTS_RETAIN_PRIME.md](RESULTS_RETAIN_PRIME.md)**

### `retain_prime_compact_method.py`

تحل مشكلة ذاكرة `set` مع الحفاظ على الفكرة نفسها:

- العدد `2` يعالج منفردًا.
- تمثل الأعداد الفردية فقط.
- كل مرشح فردي يأخذ **بتًا واحدًا**.
- يبدأ حذف مضاعفات `p` من `p²` وبخطوة `2p`.
- الناتج النهائي في هذه النسخة ما زال `list[int]` تقليدية.

عند `N = 100000` تحتاج مصفوفة المرشحين نفسها إلى `6250` بايت فقط.

```bash
python retain_prime_compact_method.py 100
```

التقرير:

**[RESULTS_COMPACT_MEMORY.md](RESULTS_COMPACT_MEMORY.md)**

### `retain_prime_packed_output_method.py`

تضغط **فضاء المرشحين والناتج معًا**:

- المرشحون: بت واحد لكل عدد فردي.
- الناتج الكامل: `array('I')` أو `array('Q')` بدل `list[int]`.
- عند عدم الحاجة إلى الاحتفاظ بكل النتائج: `iter_primes_retain_packed()` كواجهة streaming.

```bash
python retain_prime_packed_output_method.py 100
python retain_prime_packed_output_method.py 100 --stream
```

التقرير:

**[RESULTS_PACKED_OUTPUT.md](RESULTS_PACKED_OUTPUT.md)**

### `segmented_method.py`

المرحلة التالية للتوسع في `N`: بدل تخصيص bitset لكل الأعداد الفردية حتى `N`، يقسم المجال إلى **مقاطع ثابتة السعة**.

السعة الافتراضية:

```text
32768 عددًا فرديًا = 4096 bytes للمقطع
```

في كل مقطع:

- يبقى العدد الأولي ولا يُحذف.
- تحذف المركبات باستخدام أوليات الأساس حتى `sqrt(N)`.
- يبدأ شطب `p` من `p²` أو أول مضاعف له داخل المقطع.
- تكون خطوة الشطب `2p` لأن الزوجيات غير ممثلة.
- تُخرج النتائج streaming ثم يُتخلص من المقطع.

واجهات الاستخدام:

```python
iter_primes_segmented(...)
consume_primes_segmented(...)
primes_segmented_packed(...)
```

تشغيل مثال:

```bash
python segmented_method.py 1000000
python segmented_method.py 1000000 --segment-odds 32768
```

التقرير الرسمي:

**[RESULTS_SEGMENTED.md](RESULTS_SEGMENTED.md)**

### `optimized_method.py`

نسخة غربال محسنة تستخدم `bytearray`، تبدأ الشطب من `p²`، وتتجنب المرور على المضاعفات الزوجية بعد معالجة `2`.

```bash
python optimized_method.py 100
```

### `trace_method.py`

نسخة تعليمية تعرض القوائم **أ، ب، ج...** والمضاعفات المحذوفة والأعداد الأولية الجديدة مرحلةً بمرحلة.

```bash
python trace_method.py 100
```

مثال كامل:

**[TRACE_100.md](TRACE_100.md)**

### `ci_verify.py`

مرجع تحقق مستقل يعتمد القسمة التجريبية حتى `sqrt(n)`، ويقارن جميع التطبيقات، بما فيها `array` والـstreaming والتجزئة.

```bash
python ci_verify.py 4999
```

## المقارنة التجريبية

### المقارنة العامة

`benchmark.py` يقيس:

- صحة النتائج.
- زمن التنفيذ.
- ذروة الذاكرة بواسطة `tracemalloc`.
- تمثيلات الناتج المختلفة: `list` و`array` وstreaming.

```bash
python benchmark.py
python benchmark.py 1000 10000 100000 1000000 --repeats 5
python benchmark.py 1000 10000 100000 --repeats 5 --csv results.csv
```

### مقارنة التجزئة

`benchmark_segmented.py` يقارن **streaming كامل المجال** مع **streaming مجزأ** من دون الاحتفاظ بقائمة ناتج في أي من الطرفين.

```bash
python benchmark_segmented.py
python benchmark_segmented.py 100000 1000000 --repeats 5 --segment-odds 32768
```

ويتحقق من التطابق بواسطة عدد الأوليات ومجموعها وآخر عدد أولي.

## النتائج الموثقة

- **[RESULTS.md](RESULTS.md)** — خط الأساس الأولي.
- **[RESULTS_RETAIN_PRIME.md](RESULTS_RETAIN_PRIME.md)** — تجربة إبقاء العدد الأولي.
- **[RESULTS_GITHUB_ACTIONS.md](RESULTS_GITHUB_ACTIONS.md)** — أول benchmark رسمي على GitHub Actions.
- **[RESULTS_COMPACT_MEMORY.md](RESULTS_COMPACT_MEMORY.md)** — ضغط فضاء المرشحين إلى بت واحد لكل عدد فردي.
- **[RESULTS_PACKED_OUTPUT.md](RESULTS_PACKED_OUTPUT.md)** — ضغط الناتج بواسطة `array` واختبار streaming.
- **[RESULTS_SEGMENTED.md](RESULTS_SEGMENTED.md)** — تقسيم فضاء البحث إلى مقاطع ثابتة ودراسة توسع الذاكرة حتى `N = 10,000,000`.

## خلاصة ضغط الذاكرة عند N = 100000

في Benchmark رسمي على GitHub Actions:

| التمثيل | ذروة الذاكرة |
|---|---:|
| إبقاء الأولي باستخدام `set` | 8,800,704 B |
| النسخة المحسنة السابقة | 491,206 B |
| bitset + `list[int]` | 397,491 B |
| **bitset + `array`** | **47,767 B** |
| **bitset + streaming** | **7,767 B** |

عند هذا الحد، تخزين الناتج في `array` خفض الذاكرة بنحو **87.98%** مقارنة بالنسخة التي تستخدم bitset لكنها تعيد `list[int]`، وبنحو **90.28%** مقارنة بالنسخة المحسنة السابقة.

## خلاصة التجزئة والتوسع

بمقاطع سعتها `32768` عددًا فرديًا:

| N | full streaming peak | segmented streaming peak | نسبة full / segmented |
|---:|---:|---:|---:|
| 100,000 | 7,787 B | 8,506 B | 0.92× |
| 1,000,000 | 64,029 B | **10,872 B** | **5.89×** |
| 10,000,000 | 626,625 B | **12,116 B** | **51.72×** |

عند `N = 100000` لا تفيد التجزئة لأن bitset الكامل صغير أصلًا. لكن عند `N = 10000000` خفضت الذاكرة بنحو **98.07%** مقارنة بالـstream الكامل.

والأهم: عندما زاد `N` بمقدار `100×` من `100000` إلى `10000000`، زادت ذروة الذاكرة المجزأة فقط من نحو `8.5 KB` إلى `12.1 KB`.

هذه الذاكرة ليست ثابتة رياضيًا بالكامل، لأن أوليات الأساس حتى `sqrt(N)` ما زالت تنمو مع `N`. لكن ذاكرة **فضاء البحث الرئيسي** أصبحت مرتبطة بحجم المقطع بدل `N` كله.

## كلفة الزمن في النسخة المجزأة

التحكم القوي في الذاكرة له كلفة في تطبيق Python الحالي. في القياسات الرسمية كانت التجزئة أبطأ بنحو `1.7×` من streaming الكامل.

عند `N = 10,000,000`، بتكرار واحد لفحص اتجاه التوسع:

```text
full stream = 46.066240 s
segmented   = 78.815504 s
```

الغرض من النسخة المجزأة حاليًا هو **التوسع تحت قيد ذاكرة صغير**، لا تحقيق أفضل زمن ممكن.

## اختبارات الصحة المحلية

```bash
python -m unittest -v
```

تتحقق الاختبارات من:

- النتائج المعروفة حتى 100.
- الحالات الحدية الصغيرة.
- تطابق جميع الطرق حتى نطاقات اختبار متعددة.
- صحة تمثيل البتات.
- صحة `array` المضغوطة.
- تطابق streaming مع المرجع.
- عبور حدود مقاطع متعددة في الطريقة المجزأة.
- أن العدد `1` لا يدخل قائمة الأعداد الأولية.

## GitHub Actions: المرجع الرسمي للتحقق

### Tests

`.github/workflows/tests.yml` يعمل تلقائيًا عند كل `push` إلى `main` وكل Pull Request، على:

```text
Python 3.11
Python 3.12
Python 3.13
```

ويشغل الترجمة، الاختبارات الوحدوية، والتحقق المستقل حتى `N = 4999`.

### Benchmark

`.github/workflows/benchmark.yml` مسار يدوي (`workflow_dispatch`) للمقارنة العامة. يسجل بيئة التشغيل والـcommit ويحفظ:

```text
benchmark-results.csv
benchmark-output.txt
environment.txt
```

### Segmented Benchmark

`.github/workflows/segmented-benchmark.yml` مسار يدوي مستقل لدراسة توسع الذاكرة مع التجزئة. يسمح بتحديد:

- حدود `N`.
- عدد التكرارات.
- عدد المرشحين الفرديين في كل مقطع.

ويحفظ:

```text
segmented-results.csv
segmented-output.txt
segmented-environment.txt
```

كـ GitHub Actions artifacts لمدة 90 يومًا.

## هدف التجربة

الهدف ليس مجرد إنتاج قائمة الأعداد الأولية، بل دراسة أثر كل قرار في الخوارزمية وتمثيل البيانات على الصحة والزمن والذاكرة، مع قياسات قابلة لإعادة الإنتاج وتوثيق كل مرحلة بصورة مستقلة.
