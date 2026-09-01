# Normalized-u fine-spectrum transfer — SVG index

التقرير المرتبط:

- [`../NORMALIZED_U_SPECTRUM_TRANSFER_STUDY.tex`](../NORMALIZED_U_SPECTRUM_TRANSFER_STUDY.tex)

## 1. الطيف المرجعي في الإحداثي المحلي المطبع

![Normalized-u reference spectrum](normalized_u_reference_spectrum.svg)

- [`normalized_u_reference_spectrum.svg`](normalized_u_reference_spectrum.svg)

كل فترة بين مربعي أوليين متتاليين لها طول 1 في الإحداثي `u`. يبين الرسم أحد حلول NNLS للطيف المرجعي؛ القيم الفردية ليست فريدة بالضرورة، ولذلك لا تُفسر كقانون مستقل.

## 2. مقارنة النقل القديم بالنقل بعد التطبيع المحلي

![Normalized-u transfer improvement](normalized_u_transfer_improvement.svg)

- [`normalized_u_transfer_improvement.svg`](normalized_u_transfer_improvement.svg)

النقل في الإحداثي الفيزيائي باستخدام `L/H_median` كان يعطي متوسط kernel RMSE يقارب `2.60e-2`. بعد تحويل كل فترة إلى وحدة واحدة في `u` انخفض متوسط الخطأ على 17 لوحة غير مستخدمة في الملاءمة إلى نحو `4.88e-4` في النسخة uniform، أي تحسن بعامل يقارب 53.

## 3. نقل spectrum واحد من 10^7 إلى 10^12

![Normalized-u high-X kernel transfer](normalized_u_high_X_kernel_transfer.svg)

- [`normalized_u_high_X_kernel_transfer.svg`](normalized_u_high_X_kernel_transfer.svg)

الخط المتصل هو target kernel النظري في ثلاث لوحات قرب `2e12`، والمتقطع هو prediction الناتج من spectrum واحد ملائم قرب `2e7` مع تعديل amplitude اللوغاريتمي فقط ومن دون إعادة ملاءمة شكل الطيف.
