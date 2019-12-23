################
Validation Tests
################

{% for test in test_json %}

.. raw:: html

    <h2>{{ test["name"] }}</h2>
    <div>
        <table>
            <tr>
                <td>id</td>
                <td><a href="" target="_blank">{{ test["id"] }}</a></td>
            </tr>
            <tr>
                <td>uri</td>
                <td>{{ test["uri"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>name</td>
                <td>{{ test["name"] }}</td>
            </tr>
            <tr>
                <td>alias</td>
                <td>{{ test["alias"] }}</td>
            </tr>
            <tr>
                <td>author</td>
                <td>
                    {% for name in test["author"] %} {{ name["given_name"] }} {{ name["family_name"] }} {{ ", " if not loop.last }} {% endfor %}
                </td>
            </tr>
            <tr>
                <td>creation_date</td>
                <td>{{ test["creation_date"] }}</td>
            </tr>
            <tr>
                <td>status</td>
                <td><span> {{ test["status"] }}</span></td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>species</td>
                <td>{{ test["species"] }}</td>
            </tr>
            <tr>
                <td>brain_region</td>
                <td>{{ test["brain_region"] }}</td>
            </tr>
            <tr>
                <td>cell_type</td>
                <td>{{ test["cell_type"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>data_location</td>
                <td><a href='{{ test["data_location"] }}' target="_blank">{{ test["data_location"] }}</a></td>
            </tr>
            <tr>
                <td>data_type</td>
                <td>{{ test["data_type"] }}</td>
            </tr>
            <tr>
                <td>data_modality</td>
                <td>{{ test["data_modality"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>test_type</td>
                <td>{{ test["test_type"] }}</td>
            </tr>
            <tr>
                <td>score_type</td>
                <td>{{ test["score_type"] }}</td>
            </tr>
            <tr class="brown lighten-5">
                <td colspan="2"></td>
            </tr>
            <tr>
                <td>protocol</td>
                <td>{{ test["protocol"].encode('unicode_escape')|e }}</td>
            </tr>
        </table>
    </div>

{% endfor %}