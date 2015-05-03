#include "Light.h"

void Light::setColor(Vector3 v){
    Vector3 ppa = m_material->powerPerArea();
    m_material->setPowerPerArea(ppa*v / m_color);
    m_color = v;
}

void Light::setWattage(float f) {
    Vector3 ppa = m_material->powerPerArea();
    m_material->setPowerPerArea(ppa*f / m_wattage);
    m_wattage = f;
}