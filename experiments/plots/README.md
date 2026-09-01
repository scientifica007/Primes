# الرسوم الموثقة للتجارب

هذا المجلد مخصص للرسوم التي تنتج عن التجارب الموثقة في `experiments/`.

## تجربة نموذج الارتباط الزوجي

التقرير المرتبط:

- [`../PAIR_CORRELATION_MODEL_EXPERIMENT.md`](../PAIR_CORRELATION_MODEL_EXPERIMENT.md)

الرسوم:

### 1. مقارنة النموذج الزوجي الأول مع نموذج الاستقلال

![First pair-correlation model](pair_correlation_block_variance_comparison.svg)

الملف:

- [`pair_correlation_block_variance_comparison.svg`](pair_correlation_block_variance_comparison.svg)

### 2. توسيع مدى علاقات الأزواج حتى 50 فترة

![Pair-correlation range sweep](pair_correlation_range_sweep.svg)

الملف:

- [`pair_correlation_range_sweep.svg`](pair_correlation_range_sweep.svg)

### 3. توزيع معاملات الارتباط النظرية للأزواج المتجاورة

![Theoretical adjacent pair correlations](pair_correlation_rho_distribution.svg)

الملف:

- [`pair_correlation_rho_distribution.svg`](pair_correlation_rho_distribution.svg)

## تجربة نموذج Hardy–Littlewood النقطي بعد عجلة 997

التقرير المرتبط:

- [`../HL_SINGULAR_SERIES_POINT_MODEL.md`](../HL_SINGULAR_SERIES_POINT_MODEL.md)

الرسوم:

### 1. تطور الارتباط المتتالي مع تعميق الغربلة العشوائية

![HL point-model lag-1 progression](hl_point_model_lag1_progression.svg)

الملف:

- [`hl_point_model_lag1_progression.svg`](hl_point_model_lag1_progression.svg)

يبين أن ذيل `singular series` بعد عجلة 997 يدفع الارتباط في الاتجاه السلبي الصحيح، لكنه لا يصل إلى الارتباط المرصود.

### 2. تطور التشتت مع تعميق الغربلة العشوائية

![HL point-model dispersion progression](hl_point_model_dispersion_progression.svg)

الملف:

- [`hl_point_model_dispersion_progression.svg`](hl_point_model_dispersion_progression.svg)

يبين المقايضة التي ظهرت في التجربة: الارتباط يتحسن في اتجاه البيانات، بينما يصبح التشتت أصغر من الواقع مع زيادة `Q`.

### 3. فحص مباشر لوزن الأزواج المتبقي بعد عجلة 997

![Residual singular-series pair-weight check](hl_residual_pair_weight_check.svg)

الملف:

- [`hl_residual_pair_weight_check.svg`](hl_residual_pair_weight_check.svg)

يقارن النسبة التجريبية لأزواج الأوليات بين ناجي عجلة 997 بالوزن النظري المتبقي لـ`singular series` حسب المسافة `h`.

## مقارنة عجلات الأساس 29 و97 و251 و997

التقرير الرياضي محفوظ بصيغة LaTeX لتجنب مشاكل عرض التعابير المعقدة في Markdown:

- [`../BASE_WHEEL_PAIR_BALANCE_STUDY.tex`](../BASE_WHEEL_PAIR_BALANCE_STUDY.tex)

### 1. مسح مقياس التوازن عبر قيم Q

![Base-wheel balance sweep](base_wheel_pair_balance_score.svg)

- [`base_wheel_pair_balance_score.svg`](base_wheel_pair_balance_score.svg)

يبين كيف يتغير مقياس التوازن بين التشتت والارتباط عند تثبيت عجلات أساس مختلفة ثم تعميق ذيل الغربلة العشوائية.

### 2. جولات التأكيد لأفضل الإعدادات

![Confirmed base-wheel balance](base_wheel_pair_balance_confirmation.svg)

- [`base_wheel_pair_balance_confirmation.svg`](base_wheel_pair_balance_confirmation.svg)

النقطة المثالية في هذا الرسم هي الأصل `(0,0)`: صفر خطأ معياري في التشتت وصفر خطأ معياري في الارتباط. أفضل إعداد وفق مقياس المسافة المستخدم كان `q0=251, Q=7919`، مع منافسة قوية من `q0=29, Q=7919`.

## المسح الموسع لعجلة الأساس عند Q=7919

التقرير الرياضي:

- [`../BASE_WHEEL_Q7919_SWEEP_STUDY.tex`](../BASE_WHEEL_Q7919_SWEEP_STUDY.tex)

### 1. منحنى خطأ التشتت وخطأ الارتباط

![Fixed-Q discrepancy curve](base_wheel_Q7919_discrepancy_curve.svg)

- [`base_wheel_Q7919_discrepancy_curve.svg`](base_wheel_Q7919_discrepancy_curve.svg)

يعرض كيف يتغير خطأ التشتت وخطأ الارتباط عندما نزيد مقدار البنية المثبتة داخل عجلة الأساس مع تثبيت نهاية الذيل عند `Q=7919`.

### 2. منحنى مقياس التوازن عبر عجلات الأساس

![Fixed-Q balance curve](base_wheel_Q7919_balance_curve.svg)

- [`base_wheel_Q7919_balance_curve.svg`](base_wheel_Q7919_balance_curve.svg)

المسح لا يظهر عجلة سحرية واحدة؛ توجد هضبة واسعة من الإعدادات المقبولة، بينما تتدهور المطابقة بوضوح عند عجلات الأساس الكبيرة جدًا. نقاط الدائرة المفتوحة هي الإعدادات التي أُعيدت بمحاكاة أكبر للتأكيد.

## المقارنة المباشرة مع Leung 2026 وHardy–Littlewood

التقرير الرياضي:

- [`../THEORY_ALIGNED_LEUNG_HL_COMPARISON.tex`](../THEORY_ALIGNED_LEUNG_HL_COMPARISON.tex)

هذه المرحلة تختبر أولًا ما إذا كانت `pair covariance` المعروفة في نظرية الفترات القصيرة تفسر السلوك الذي كنا على وشك نسبته إلى علاقات ثلاثية مستقلة.

### 1. الارتباط في العد الموزون مقارنة بتنبؤ Leung 2026

![Leung correlation comparison](leung2026_correlation_comparison.svg)

- [`leung2026_correlation_comparison.svg`](leung2026_correlation_comparison.svg)

يعرض معاملات الارتباط المرصودة في أربع فترات قصيرة متجاورة ذات طول نسبي ثابت، مقارنة بالتنبؤ المستخرج من صيغة التغاير في Leung 2026. عند القيم الأصغر لـ`delta` تتفق الإشارة والحجم جيدًا، ويقلل الحد الثانوي جزءًا مهمًا من انحياز الحجم المنتهي.

### 2. سلوك runs المحلي: البيانات مقابل Gaussian تحدده pair covariance

![Leung runs comparison](leung2026_runs_comparison.svg)

- [`leung2026_runs_comparison.svg`](leung2026_runs_comparison.svg)

يبين أن متوسط عدد تبدلات الإشارة في أربع فترات متجاورة قريب من نموذج Gaussian متعدد المتغيرات الذي تحدده مصفوفة الأزواج وحدها.

### 3. أنماط الإشارات الثلاثية في فترات مربعات الأوليات

![Prime-square triple sign patterns](prime_square_pair_gaussian_sign_triples.svg)

- [`prime_square_pair_gaussian_sign_triples.svg`](prime_square_pair_gaussian_sign_triples.svg)

يقارن الأنماط الثمانية `000` إلى `111` في 448 ثلاثية متجاورة مع نموذج pair-Gaussian. لم يظهر أي نمط ثلاثي شاذ إحصائيًا، وهو سبب مباشر لتأجيل الانتقال إلى نموذج `k=3` صريح.

## النموذج النقطي لمنحنى الارتباط الكامل rho(k)

التقرير الرياضي:

- [`../POINT_MODEL_FULL_RHO_KERNEL_STUDY.tex`](../POINT_MODEL_FULL_RHO_KERNEL_STUDY.tex)

هذه التجربة تشخّص أولًا أن Bernoulli thinning المستقل يضيف تباينًا بلا تغاير فيضعف الارتباط، ثم تستبدله بتجربة block-quota، وأخيرًا تضيف حقلًا زوجيًا واحدًا ذا covariance محددة من نظرية الفترات القصيرة. لم تُستخدم أي معلومة ثلاثية في الملاءمة.

### 1. منحنى الارتباط الكامل حتى lag 20

![Full rho kernel](point_model_full_rho_kernel_comparison.svg)

- [`point_model_full_rho_kernel_comparison.svg`](point_model_full_rho_kernel_comparison.svg)

يبين أن random sieve الخام يحمل بالفعل kernel زوجيًا قريبًا من النظرية، وأن Pair-kernel point model يحافظ على هذا الشكل مع تصحيح التباين الهامشي.

### 2. مقارنة شكل kernel مباشرة بدالة Delta(k)

![Delta kernel shape](point_model_delta_shape_comparison.svg)

- [`point_model_delta_shape_comparison.svg`](point_model_delta_shape_comparison.svg)

يعرض `rho(k)/rho(1)` مقابل `Delta(k)/Delta(1)`. خطأ الشكل المطبع للنموذج النهائي يقارب `0.037`.

### 3. المقايضة بين التشتت والارتباط

![Variance-correlation tradeoff](point_model_variance_correlation_tradeoff.svg)

- [`point_model_variance_correlation_tradeoff.svg`](point_model_variance_correlation_tradeoff.svg)

يوضح لماذا لا يمكن لضوضاء thinning مستقلة أن تطابق التشتت والارتباط معًا: رفع التباين بهذه الطريقة يخفض مقدار الارتباط. الحقل الزوجي المصحح ينتقل إلى قرب التشتت المرصود من دون التضحية بالـkernel.

## الآلية النقطية المحلية متعددة المقاييس

التقرير الرياضي:

- [`../LOCAL_MULTISCALE_PAIR_MECHANISM_STUDY.tex`](../LOCAL_MULTISCALE_PAIR_MECHANISM_STUDY.tex)

هذه المرحلة تنقل تصحيح `pair covariance` من حقل Gaussian على مستوى الفترات إلى عمليات إعادة توزيع محلية محدودة الدعم على خط الأعداد، باستخدام موجات Haar متعددة المقاييس. أوزان المقاييس اختيرت من الـkernel النظري، بينما استُخدم من بيانات التدريب معامل شدة واحد فقط لضبط مقدار التشتت.

### 1. أوزان المقاييس المحلية

![Local multiscale weights](local_multiscale_scale_weights.svg)

- [`local_multiscale_scale_weights.svg`](local_multiscale_scale_weights.svg)

يوضح المقاييس المحلية التي اختارها NNLS لتمثيل الـkernel النظري. ظهر أكبر وزن عند دعم يساوي `65536` عددًا صحيحًا في هذه العينة، لكن لا ينبغي تفسير هذا الرقم كحد عالمي أو ثابت قبل اختبار قانون scaling عبر قيم مختلفة لـ`X` و`H`.

### 2. منحنى rho(k) الكامل

![Local multiscale full rho](local_multiscale_full_rho_comparison.svg)

- [`local_multiscale_full_rho_comparison.svg`](local_multiscale_full_rho_comparison.svg)

يقارن البيانات والنظرية والنموذج النقطي المحلي عبر `lag=1..20`. النموذج المحلي يطابق التشتت وقوة `rho(1)` جيدًا، والمنحنى المرصود كاملًا لا يبدو شاذًا تحت المحاكاة.

### 3. مقارنة الشكل بدالة Delta(k)

![Local multiscale Delta shape](local_multiscale_delta_shape_comparison.svg)

- [`local_multiscale_delta_shape_comparison.svg`](local_multiscale_delta_shape_comparison.svg)

يبين أن الآلية المحلية تولد amplitude قريبًا جدًا من النظرية وشكلًا قريبًا من `Delta(k)` من دون إدخال مصفوفة covariance Gaussian جاهزة على مستوى الفترات. بقي خطأ الشكل أكبر من bridge model، ولذلك الخطوة التالية هي البحث عن قانون scaling لأوزان المقاييس بدل إعادة ملاءمتها عدديًا في كل نطاق.

## قاعدة التوثيق

عندما تنتج تجربة رسمًا يستخدم في الاستنتاج أو المقارنة، يجب حفظ نسخة قابلة للعرض داخل المستودع وربطها بالتقرير أو بهذا الفهرس، حتى لا تبقى الرسوم محصورة في ملفات الجلسة المحلية.
