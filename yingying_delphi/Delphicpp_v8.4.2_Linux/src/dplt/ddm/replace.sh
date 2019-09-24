find . -type f -exec sed -ri '/(^#include <ddm)/s/</\"..\//' {} \; -exec sed -ri '/(^#include ")/s/>/\"/' {} \; -exec sed -ri '/(^#include ")/s/ddm\/dart\/if/dart-impl/' {} \;
