#include "Light.h"

void Light::setColor(Vector3 v){
    m_color = v;
    m_material->setPowerPerArea(v);
}