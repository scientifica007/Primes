# خريطة الأعمال السابقة وموقع المشروع

تاريخ المراجعة: **31 أغسطس 2026**.

## لماذا هذه الوثيقة؟

هدفها منع المشروع من الخلط بين ثلاثة أشياء مختلفة:

1. **إعادة اكتشاف أو إعادة تنفيذ فكرة معروفة**.
2. **نتيجة هندسية خاصة بتطبيقنا وبيئة القياس**.
3. **سؤال بحثي قد يبقى مفتوحًا ويستحق متابعة مستقلة**.

من الآن فصاعدًا لا تُوصف أي نتيجة بأنها جديدة قبل مقارنتها بهذه الخريطة وتحديثها إن لزم.

---

## الخلاصة التنفيذية

معظم المسار الخوارزمي الذي وصلنا إليه معروف في أدبيات الغربلة:

- الغربلة المجزأة Segmented Sieve معروفة على الأقل منذ Bays وHudson (1977).
- فكرة تجنب مضاعفات الأوليات الصغيرة بعجلة Wheel صيغت ودرست بوضوح عند Paul Pritchard في 1981–1983.
- المقايضات بين الزمن والذاكرة في غرابيل الأوليات درست عند Jonathan Sorenson، خصوصًا في 1990 و1998.
- تصميم غربلة مجزأة cache-friendly مع bucket scheduling للأوليات الكبيرة طوّره Tomás Oliveira e Silva بدءًا من 2001 ووثقه في 2002/2003.
- `primesieve` الحديث يطبق segmented Eratosthenes + wheel factorization + cache-size tuning + bucket sieve وتحسينات SIMD/branch/cache إضافية.
- Atkin–Bernstein (2004) يمثل مسارًا نظريًا مختلفًا مبنيًا على binary quadratic forms.

إذن **مشروعنا حتى الآن ليس خوارزمية أوليات جديدة**. قيمته الحالية هي: إعادة بناء منهجية، قياسات قابلة لإعادة الإنتاج، ودراسة تجريبية لتحديد أين توجد أسئلة لا تختزل إلى prior art.

---

## جدول المطابقة

| فكرة في مشروعنا | أقرب prior art | الحالة |
|---|---|---|
| حذف مضاعفات الأوليات / البدء من `p²` | Sieve of Eratosthenes وتحسيناته القياسية | **KNOWN** |
| عدم تمثيل الزوجيات | تحسين كلاسيكي | **KNOWN** |
| bitset بدل كائن لكل عدد | تطبيقات الغربلة المدمجة الحديثة | **KNOWN** |
| segmented sieve | Bays & Hudson, 1977 | **KNOWN** |
| تمثيل فقط الأعداد المتباينة أوليًا مع primorial | Pritchard wheel sieve, 1981–1983 | **KNOWN** |
| `6 -> 30 -> 210 -> 2310 -> 30030` | wheel/primorial sieves | **KNOWN** |
| الكثافة `phi(M)/M` | Euler totient / wheel theory | **KNOWN** |
| مقايضة wheel أكبر مقابل جدول أكبر | wheel-sieve engineering literature | **KNOWN** |
| اختيار حجم المقطع وفق L1/L2 | Oliveira e Silva; primesieve | **KNOWN** |
| bucket scheduling للأوليات الكبيرة | Oliveira e Silva, 2001+ | **KNOWN** ولم نطبقه بعد |
| جداول wheel مضغوطة fixed-width | تحسين تمثيل هندسي طبيعي؛ توجد أشكال عديدة منه | **REPLICATION/ENGINEERING** |
| نتائج أفضل حجم مقطع على GitHub runner | خاص بتطبيقنا/المعالج/بايثون | **PROJECT-SPECIFIC BENCHMARK** |
| فحص آخر رقم/آخر رقمين/ثلاثة في base `b` | arithmetic progressions modulo `b^k` | **KNOWN FRAMEWORK** |
| مجموع الأرقام وعلاقته بـ `b-1` | اختبارات القسمة + أدبيات sum-of-digits | **KNOWN** |
| المجموع المتناوب وعلاقته بـ `b+1` | اختبار divisibility modular تقليدي | **KNOWN** |
| equidistribution لمجموع أرقام الأوليات | Mauduit–Rivat, 2010 | **KNOWN DEEP RESULT** |
| انتقالات بواقي **الأوليات المتتالية** | Lemke Oliver–Soundararajan, 2016 | **KNOWN NONTRIVIAL BIAS** |
| انتقالات الأرقام داخل **العدد الأولي نفسه** بعد إزالة modular effects | لم نجد في هذا المسح نتيجة مطابقة مباشرة | **OPEN/NEEDS LITERATURE SEARCH** |
| مقياس معلوماتي يبحث عن base يميز prime/composite بعد conditioning الكامل على congruences | لم نجد نتيجة مطابقة مباشرة | **OPEN/EXPERIMENTAL** |
| اختيار wheel تكيفية آليًا من نموذج cache/memory لا من benchmark فقط | توجد auto-tuning أفكار مشابهة؛ الجدة غير مثبتة | **OPEN ONLY IF FORMALIZED** |

---

## 1. Segmented Sieve

Carter Bays وRichard H. Hudson نشرا في 1977:

> *The segmented sieve of Eratosthenes and primes in arithmetic progressions to 10^12*, BIT 17, 121–127.

DOI: https://doi.org/10.1007/BF01932283

هذا يضع فكرة تقسيم المجال إلى segments خارج أي ادعاء جدة في مشروعنا.

---

## 2. Pritchard والـWheel

أعمال Paul Pritchard الأساسية:

1. **1981** — *A Sublinear Additive Sieve for Finding Prime Numbers*, Communications of the ACM 24(1), 18–23.  
   DOI: https://doi.org/10.1145/358527.358540
2. **1982** — *Explaining the Wheel Sieve*, Acta Informatica 17(4), 477–485.  
   DOI: https://doi.org/10.1007/BF00264164
3. **1983** — *Fast Compact Prime Number Sieves (among Others)*, Journal of Algorithms 4(4), 332–344.  
   DOI: https://doi.org/10.1016/0196-6774(83)90014-7

ورقة 1983 تصف عائلة parameterized من الغرابيل وتتعامل صراحةً مع الاقتصار على الأعداد غير القابلة للقسمة على أول `k` أوليات. هذه هي البنية الرياضية الأساسية لما أسميناه في المشروع:

```text
Wheel-6
Wheel-30
Wheel-210
Wheel-2310
Wheel-30030
```

إذن تسلسل primorial نفسه ليس جديدًا.

---

## 3. Sorenson: تحليل الغربلة ومقايضة الزمن/الذاكرة

Jonathan Sorenson:

- *An Introduction to Prime Number Sieves*, University of Wisconsin Technical Report TR909, 1990.  
  https://minds.wisconsin.edu/handle/1793/59248
- *Trading Time for Space in Prime Number Sieves*, ANTS-III, 1998, pp. 179–195.  
  DOI: https://doi.org/10.1007/BFb0054861

هذه الأعمال مهمة لأن مشروعنا مر مرارًا بنفس السؤال: هل نزيد wheel/الضغط/التجزئة مقابل زمن إضافي؟ هذا النوع من space-time tradeoff جزء صريح من الأدبيات، وليس اتجاهًا جديدًا بحد ذاته.

---

## 4. Oliveira e Silva: cache وbucket scheduling

Tomás Oliveira e Silva وثق *Fast implementation of the segmented sieve of Eratosthenes*، وهي خوارزمية صممت أصلًا لحسابات قرب `10^17` و`10^18`.

المصدر: https://sweet.ua.pt/tos/software/prime_sieve.html

النقاط ذات الصلة المباشرة بمشروعنا:

- bit واحد لكل مرشح بعد pre-sieving بعوامل صغيرة مثل `2,3,5`.
- segments صغيرة لتبقى البيانات المهمة في cache.
- مشكلة أن المرور على **كل sieving primes لكل segment** يصبح مكلفًا عندما يكون `N` ضخمًا.
- الحل: توزيع sieving primes على قوائم/buckets مرتبطة بالمقاطع التي ستضرب فيها لاحقًا، ثم نقلها إلى bucket التالي عند المعالجة.

هذه نقطة يتوقف عندها تطبيقنا الحالي: نحن ما زلنا نمر على base primes في كل segment ثم نقرر متى نتوقف. إذا أردنا تجاوز prior art الحالي فعلينا أولًا تنفيذ bucket scheduling المعروف ومقارنته، لا الادعاء أن تقليل ذاكرة segment وحده جديد.

---

## 5. primesieve كـbaseline احترافي

المشروع: https://github.com/kimwalisch/primesieve

وثيقة الخوارزميات: https://github.com/kimwalisch/primesieve/blob/master/doc/ALGORITHMS.md

`primesieve` يستخدم:

- segmented Sieve of Eratosthenes؛
- wheel factorization؛
- bit array؛
- اكتشاف L1/L2 واختيار sieve size بناءً عليه؛
- خوارزميات مختلفة للأوليات الصغيرة والمتوسطة والكبيرة؛
- bucket sieve مستمد من Oliveira e Silva عند الحاجة؛
- SIMD وتحسينات لتقليل branch misprediction واستخدام cache hierarchy.

من الآن فصاعدًا **السرعة مقارنة بإصدار أقدم من كودنا لا تكفي**. أي ادعاء أداء يجب أن يتضمن baseline خارجيًا مثل `primesieve`، ويفضل single-thread عند مقارنة خوارزمية single-thread.

---

## 6. Atkin–Bernstein: مسار مختلف نظريًا

A. O. L. Atkin وDaniel J. Bernstein:

*Prime sieves using binary quadratic forms*, Mathematics of Computation 73 (2004), 1023–1030.  
DOI: https://doi.org/10.1090/S0025-5718-03-01501-1

الخوارزمية لا تكتفي بإعادة ترتيب Eratosthenes/Wheel بل تعتمد تمثيلات بأشكال تربيعية ثنائية. الورقة تعطي حدودًا نظرية مختلفة للعمليات والذاكرة.

هذا مرجع مهم لتذكيرنا بأن تحسين Wheel داخل Eratosthenes ليس المجال الوحيد لخوارزميات توليد الأوليات.

---

## 7. تغيير قاعدة العد: ما المعروف؟

### 7.1 آخر `k` أرقام

آخر `k` أرقام للعدد `n` في base `b` تحدد `n mod b^k`. لذلك توزيع suffixes للأوليات يدخل مباشرة في نظرية الأوليات في المتتاليات الحسابية/reduced residue classes.

هذا يعني أن مجرد ملاحظة أن بعض suffixes لا تظهر أو أن المسموح منها يقترب من الانتظام لا يمثل بنية رقمية جديدة.

### 7.2 مجموع الأرقام

Christian Mauduit وJoël Rivat أثبتا في 2010 نتيجة عميقة عن مجموع أرقام الأعداد الأولية في قاعدة `q >= 2`:

*Sur un problème de Gelfond : la somme des chiffres des nombres premiers*, Annals of Mathematics 171 (2010), 1591–1646.  
DOI: https://doi.org/10.4007/annals.2010.171.1591

بصورة مختصرة: مجموع أرقام الأوليات في قاعدة ثابتة يتوزع توزيعًا منتظمًا في المتتاليات الحسابية، باستثناء حالات degeneracy المعروفة.

إذن قسم digit-sum في تجربتنا يجب أن يبدأ من هذه النتيجة، لا من الصفر.

### 7.3 concatenated prime digits

Copeland وErdős أثبتا سنة 1946 أن العدد الناتج من concatenation للأوليات هو normal في القاعدة المناسبة تحت شروط عامة.

A. H. Copeland & P. Erdős, *Note on Normal Numbers*, Bulletin of the AMS 52 (1946), 857–860.

هذا لا يجيب سؤالنا نفسه عن أرقام **كل عدد أولي منفردًا**، لكنه يثبت أن دراسة “العشوائية الرقمية للأوليات” لها تاريخ أقدم بكثير من مشروعنا.

### 7.4 انتقالات بواقي الأوليات المتتالية

Robert J. Lemke Oliver وKannan Soundararajan وجدا انحيازات غير متوقعة في أنماط بواقي **الأعداد الأولية المتتالية** modulo `q`:

*Unexpected biases in the distribution of consecutive primes*, PNAS 113(31), 2016.  
DOI: https://doi.org/10.1073/pnas.1605366113  
arXiv: https://arxiv.org/abs/1603.03720

هذا مهم جدًا عند تفسير transition matrices: توجد بالفعل biases غير تافهة في انتقالات `p_n mod q -> p_{n+1} mod q`. لكنها مسألة مختلفة عن transition بين digits متجاورة داخل تمثيل `p_n` نفسه.

---

## 8. ما الذي قد يبقى مشروعًا بحثيًا حقيقيًا؟

### A. سؤال رقمي مشروط بالكامل على modular structure

صياغة ممكنة:

> بعد conditioning على كل المعلومات التي يفرضها `n mod b^k`، و`b-1`، و`b+1`، وحجم العدد/عدد خاناته، هل يبقى statistic رقمي ذو قدرة تمييز مستقرة بين prime وcomposite تتجاوز null model مناسبًا؟

لكي يكون هذا سؤالًا علميًا يجب:

1. تعريف null model قبل رؤية النتائج.
2. تصحيح multiple testing عبر القواعد والأنماط.
3. فصل train range عن validation range.
4. اختبار الإشارة على مراتب أكبر (`10^8`, `10^9`, ...).
5. البحث عن تفسير نظري؛ لا يكفي classifier أفضل قليلًا.

الحالة: **OPEN/EXPERIMENTAL**، ولا يوجد ادعاء جدة حاليًا.

### B. نموذج تكلفة wheel/segment قابل للتنبؤ

بدل brute-force benchmarking فقط، يمكن محاولة بناء نموذج يتنبأ بأفضل:

- primorial wheel؛
- segment size؛
- table representation؛
- scheduling strategy؛

من خصائص cache والمعالج وكثافة المرشحين، ثم اختباره على معماريات متعددة.

لكن توجد auto-tuning وcache-aware implementations سابقة؛ الجدة لن تثبت إلا إذا كان النموذج مختلفًا ويفسر/يتنبأ أفضل من heuristics الموجودة.

الحالة: **POSSIBLE ENGINEERING RESEARCH**, لا ادعاء جدة بعد.

### C. ذاكرة base primes / bucket hierarchy

فكرتنا السابقة عن “تجزئة هرمية لأوليات الأساس” يجب مقارنتها أولًا بأعمال space-efficient sieves وbucket sieves. لا نعتبرها جديدة قبل إتمام هذا المسح.

الحالة: **PRIOR ART INCOMPLETE**.

---

## 9. قواعد العمل البحثية الجديدة

من هذا commit فصاعدًا:

1. كل تجربة جديدة تبدأ بقسم `Prior art checked`.
2. نستخدم baseline خارجيًا، لا نسخ مشروعنا فقط.
3. نفرق بين:
   - mathematical novelty؛
   - implementation novelty؛
   - machine-specific tuning؛
   - replication.
4. نتيجة سلبية موثقة أفضل من ادعاء إيجابي بلا null model.
5. أي “نمط” يختفي بعد conditioning modular يصنف فورًا **explained**.
6. أي benchmark لـPython أمام C/C++ لا يستخدم وحده للحكم على الخوارزمية؛ نفصل أثر اللغة عن أثر الخوارزمية.
7. لا نستخدم كلمة “اكتشاف” في README إلا مع مرجع يثبت أن السؤال لم يكن معروفًا أو مع نتيجة نظرية جديدة.

---

## 10. الخطوة التالية المعيارية

أضفنا مسار benchmark خارجيًا يقارن تطبيق المشروع مع `primesieve` على نفس GitHub runner. المقارنة تسجل نتيجتين مختلفتين عمدًا:

1. **prime counting** في `primesieve` — أقصى throughput لمهمة العد، وليس نفس workload مشروعنا.
2. **prime iteration + count + checksum + last prime** عبر C++ API لـ`primesieve` — workload أقرب إلى دالة `consume_primes_compact_spec()` في مشروعنا.

هذا يمنع مقارنة غير عادلة بين “عد فقط” و“توليد كل الأوليات ومعالجتها في Python”.

---

## مراجع أساسية

- Bays, C.; Hudson, R. H. (1977), *The segmented sieve of Eratosthenes and primes in arithmetic progressions to 10^12*, BIT 17, 121–127. https://doi.org/10.1007/BF01932283
- Pritchard, P. (1981), *A Sublinear Additive Sieve for Finding Prime Numbers*, CACM 24(1), 18–23. https://doi.org/10.1145/358527.358540
- Pritchard, P. (1982), *Explaining the Wheel Sieve*, Acta Informatica 17, 477–485. https://doi.org/10.1007/BF00264164
- Pritchard, P. (1983), *Fast Compact Prime Number Sieves (among Others)*, Journal of Algorithms 4, 332–344. https://doi.org/10.1016/0196-6774(83)90014-7
- Sorenson, J. (1990), *An Introduction to Prime Number Sieves*, TR909. https://minds.wisconsin.edu/handle/1793/59248
- Sorenson, J. (1998), *Trading Time for Space in Prime Number Sieves*. https://doi.org/10.1007/BFb0054861
- Oliveira e Silva, T., *Fast implementation of the segmented sieve of Eratosthenes*. https://sweet.ua.pt/tos/software/prime_sieve.html
- Atkin, A. O. L.; Bernstein, D. J. (2004), *Prime sieves using binary quadratic forms*. https://doi.org/10.1090/S0025-5718-03-01501-1
- Mauduit, C.; Rivat, J. (2010), *Sur un problème de Gelfond : la somme des chiffres des nombres premiers*. https://doi.org/10.4007/annals.2010.171.1591
- Lemke Oliver, R. J.; Soundararajan, K. (2016), *Unexpected biases in the distribution of consecutive primes*. https://doi.org/10.1073/pnas.1605366113
- `primesieve` algorithm documentation: https://github.com/kimwalisch/primesieve/blob/master/doc/ALGORITHMS.md
